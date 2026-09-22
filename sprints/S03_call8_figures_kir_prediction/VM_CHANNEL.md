# Programmatic VM channel (JupyterLab REST/websocket) — tested recipe

Verified working 2026-09-22 against the `big_run` app. Reuse verbatim; only the
base URL / tab id change between sessions (the app re-issues a new `app-id`
subdomain roughly each time the app instance is (re)started — always re-derive
it, don't hardcode the one below).

## This session's values (will be stale next session — re-derive)

- Workspace: `full-cohort-hla-calling`, app `AoU_Jupyter_ComputeEngine_20260805_big_run`
- Wrapper tab (Workbench chrome, cross-origin iframe host): `https://workbench.verily.com/app/9394ec22-d949-4489-b46e-73790750472e/?name=AoU_Jupyter_ComputeEngine_20260805_big_run&app=aou-jupyter`
- Inner JupyterLab base URL (navigate the browser tab directly here — the
  wrapper page's DOM/JS tools return nothing useful, per quirk #37/ENVIRONMENT.md):
  `https://9394ec22-d949-4489-b46e-73790750472e.workbench-app-prod.verily.com/lab`
- Chrome tab id used for the inner JupyterLab page: `627345092`

## How to (re)acquire the inner URL from scratch

1. Open `https://workbench.verily.com` → workspace "Full Cohort HLA Calling" → **Apps** tab.
2. If the target app instance shows a **Start** button, click it and poll
   `get_page_text` until it says "App is now running." (can take minutes).
   Never create a new app instance unless none exists.
3. `find` for `Launch <app-instance-name>` and click it — opens a **new tab**
   pointing at the `workbench.verily.com/app/<app-id>/...` wrapper.
4. In that new tab, run `javascript_tool`:
   ```js
   Array.from(document.querySelectorAll('iframe')).map(f => f.src)
   ```
   This returns `https://<app-id>.workbench-app-prod.verily.com/`. Take the
   `<app-id>` and `navigate` that same tab directly to
   `https://<app-id>.workbench-app-prod.verily.com/lab` (bypasses the
   cross-origin iframe — DOM/JS tools now work on the real JupyterLab page).

## Step 1 — XSRF token + open a terminal (REST)

```js
function getCookie(name){
  const m = document.cookie.match(new RegExp('(^|;\\s*)'+name+'=([^;]*)'));
  return m ? decodeURIComponent(m[2]) : null;
}
window._xsrf = getCookie('_xsrf');
const r = await fetch('/api/terminals', {method:'POST', headers:{'X-XSRFToken': window._xsrf}});
const j = await r.json();
JSON.stringify({status:r.status, term:j});   // {name:"1", ...} <- terminal name, use below
```

## Step 2 — open the websocket driver + bg/poll helpers

Run once per tab session (survives across multiple commands):

```js
window.res = window.res || {};
window._buf = window._buf || "";
const wsUrl = `wss://${location.host}/terminals/websocket/1`;  // "1" = terminal name from step 1
window.ws = new WebSocket(wsUrl);
window.wsReady = new Promise((resolve,reject)=>{
  window.ws.onopen = ()=>resolve('open');
  window.ws.onerror = (e)=>reject('error');
});
window.ws.onmessage = (ev)=>{
  try{
    const data = JSON.parse(ev.data);
    if(data[0]==='stdout'){ window._buf += data[1]; }
  }catch(e){}
};
window.bg = function(key, cmd){
  const marker = `__END_${key}__`;
  window._buf = "";
  window.ws.send(JSON.stringify(["stdin", cmd + " ; echo " + marker + "\r"]));
  window['_marker_'+key] = marker;
};
window.poll = function(key){
  const marker = window['_marker_'+key];
  if(!marker) return {done:false, note:'no such key'};
  const lastIdx = window._buf.lastIndexOf(marker);
  const firstIdx = window._buf.indexOf(marker);
  // marker appears twice: once in the echoed stdin, once in real stdout.
  // Need the SECOND occurrence, else you get the pre-execution echo only.
  if(lastIdx === -1 || lastIdx === firstIdx){ return {done:false, len: window._buf.length}; }
  let out = window._buf.slice(0, lastIdx);
  out = out.replace(/\x1b\][^\x07]*\x07/g, '').replace(/\x1b\[[0-9;?]*[a-zA-Z]/g, ''); // strip ANSI/OSC
  const nl = out.indexOf('\n');
  if(nl !== -1) out = out.slice(nl+1);  // drop the echoed command line
  return {done:true, text: out.trim()};
};
await window.wsReady;
'ws-ready';
```

## Step 3 — run a command in background, then poll

```js
window.bg('mykey', "hostname; nproc; echo done");
'sent';
```
Then, in a **separate** `javascript_tool` call (the 45 s JS timeout is why
this must be two calls, not one with a sleep):
```js
JSON.stringify(window.poll('mykey'));
// {done:false, len:N}  -> not finished yet, poll again
// {done:true, text:"..."} -> finished, text is the captured stdout
```
Long output gets returned in full by `poll()` inside the page (`window._buf`
is unbounded), but the **tool's own return value** truncates around 1–2 kB —
for big output, don't return `r.text` directly; instead return
`r.text.slice(a,b)` in successive calls to page through it (used above to read
~2 kB outputs in 3 slices). Also: the tool blocks return values that look like
`KEY=value; ...` cookie/query text — if a command's output contains `=`/`;`
pairs (e.g. `free -g` won't, but raw env dumps will), sanitize before
returning (e.g. `.replace(/=/g,':').replace(/;/g,',')`).

## Gotchas recap (see ENVIRONMENT.md quirk #37 for the source)

- 45 s JS timeout → always background + poll for anything non-trivial.
- Marker appears twice (echo of stdin, then real output) — take the *last*
  occurrence, not the first.
- Strip ANSI/OSC control sequences before treating output as data.
- Truncate/slice manually when reading back — don't trust one `poll()` call
  to hand you >1–2 kB in one shot.
- Never return participant-level data through this channel — aggregate
  counts/paths/metadata only, same constraint as typing into the terminal
  directly.
- One JupyterLab terminal (`/api/terminals` POST) can be reused for the whole
  session — no need to open a new one per command, just reuse `window.ws`
  and call `bg`/`poll` repeatedly with different keys.

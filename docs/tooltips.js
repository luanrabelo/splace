// ── SPLACE Tooltip Configuration ──────────────────────────────────────────────
// Centralised Tippy.js v6 theme and initialisation helpers.
// Load this file after tippy.js and before splace.js.

// ── 1. Inject custom "splace" theme CSS ───────────────────────────────────────
(function injectTheme() {
    const style = document.createElement("style");
    style.textContent = `
        /* splace theme — primary colour #323795 (splace-blue-600) */
        .tippy-box[data-theme~='splace'] {
            background-color: #323795;
            color: #fff;
            font-size: 0.82rem;
            line-height: 1.5;
            border-radius: 6px;
            box-shadow: 0 4px 16px rgba(50, 55, 149, 0.35);
        }
        .tippy-box[data-theme~='splace'] .tippy-content {
            padding: 7px 12px;
        }

        /* Arrow colours per placement */
        .tippy-box[data-theme~='splace'][data-placement^='top']    > .tippy-arrow::before { border-top-color:    #323795; }
        .tippy-box[data-theme~='splace'][data-placement^='bottom'] > .tippy-arrow::before { border-bottom-color: #323795; }
        .tippy-box[data-theme~='splace'][data-placement^='left']   > .tippy-arrow::before { border-left-color:   #323795; }
        .tippy-box[data-theme~='splace'][data-placement^='right']  > .tippy-arrow::before { border-right-color:  #323795; }

        /* HTML formatting helpers inside tooltips */
        .tippy-box[data-theme~='splace'] em       { color: #a8afe5; font-style: italic; }
        .tippy-box[data-theme~='splace'] strong   { color: #d5d9f3; font-weight: 600; }
        .tippy-box[data-theme~='splace'] code     { background: rgba(255,255,255,.15); border-radius: 3px; padding: 0 4px; font-size: .8em; }
        .tippy-box[data-theme~='splace'] a        { color: #a8afe5; text-decoration: underline; }
        .tippy-box[data-theme~='splace'] hr       { border-color: rgba(255,255,255,.2); margin: 5px 0; }
        .tippy-box[data-theme~='splace'] .tip-label { color: #a8afe5; font-size: .75rem; text-transform: uppercase; letter-spacing: .04em; }
    `;
    document.head.appendChild(style);
})();

// ── 2. Shared default options ─────────────────────────────────────────────────
const TIPPY_DEFAULTS = {
    theme: "splace",
    placement: "top",
    arrow: true,
    allowHTML: true,
    maxWidth: 340,
    duration: [150, 100],
};

// ── 3. Public helpers ─────────────────────────────────────────────────────────

/**
 * (Re-)initialise all elements that carry a `data-tippy-content` attribute.
 * Safe to call after dynamic re-renders — destroys stale instances first.
 */
window.initTooltips = function () {
    if (typeof tippy === "undefined") return;
    document.querySelectorAll("[data-tippy-content]").forEach(el => {
        if (el._tippy) el._tippy.destroy();
    });
    tippy("[data-tippy-content]", { ...TIPPY_DEFAULTS });
};

/**
 * Attach (or replace) a tooltip on a single element programmatically.
 * Used for dynamically generated elements such as gene chips and record buttons.
 *
 * @param {Element}  el       Target DOM element
 * @param {string}   content  HTML string shown inside the tooltip
 * @param {object}   [opts]   Optional overrides for TIPPY_DEFAULTS
 * @returns {object|undefined} The Tippy instance, or undefined if tippy is unavailable
 */
window.createTippy = function (el, content, opts = {}) {
    if (typeof tippy === "undefined" || !el) return;
    if (el._tippy) el._tippy.destroy();
    return tippy(el, { ...TIPPY_DEFAULTS, content, ...opts });
};

// ── 4. Auto-initialise on first load ─────────────────────────────────────────
if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", window.initTooltips);
} else {
    window.initTooltips();
}

// ==UserScript==
// @name         JP Prompt/Editor Toggle with Collapse‑All & Centered Index
// @author       https://duckduckgo.com/ running GPT-OSS 120B
// @version      1.2
// @description  Collapses all editors by default, adds an “hidden, click to show” label, and centers the input tag.
// ==/UserScript==

(function () {
    'use strict';

    const PROMPT_SELECTOR = '.jp-InputPrompt.jp-InputArea-prompt';
    const EDITOR_SELECTOR = '.jp-CodeMirrorEditor.jp-Editor.jp-InputArea-editor';
    const AFTER_LABEL = 'hidden, click to show';

    // Inject CSS for centering, the ::after label, and the collapsed state
    const style = document.createElement('style');
    style.textContent = `
        /* Center the prompt text */
        ${PROMPT_SELECTOR} {
            text-align: center;
        }

        /* Base ::after (empty) – will show only when collapsed */
        ${PROMPT_SELECTOR}::after {
            content: "";
            margin-left: 8px;
            font-size: 0.9em;
            color: #777;
            cursor: pointer;
        }

        /* Show label when the prompt has the .collapsed class */
        ${PROMPT_SELECTOR}.collapsed::after {
            content: "${AFTER_LABEL}";
        }
    `;
    document.head.appendChild(style);

    // Collapse an editor and mark its prompt
    function collapse(prompt, editor) {
        editor.style.display = 'none';
        prompt.classList.add('collapsed');
    }

    // Expand an editor and clear the prompt marker
    function expand(prompt, editor) {
        editor.style.display = 'table-cell';
        prompt.classList.remove('collapsed');
    }

    // Toggle between collapsed / expanded
    function toggle(prompt, editor) {
        if (prompt.classList.contains('collapsed')) {
            expand(prompt, editor);
        } else {
            collapse(prompt, editor);
        }
    }

    // Initialise: find all prompts, collapse their editors, and bind click handler
    function init() {
        const prompts = document.querySelectorAll(PROMPT_SELECTOR);
        prompts.forEach(prompt => {
            const editor = prompt.nextElementSibling;
            if (editor && editor.matches(EDITOR_SELECTOR)) {
                collapse(prompt, editor);                 // start collapsed
                prompt.addEventListener('click', () => toggle(prompt, editor));
            }
        });
    }

    // Run after DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }
})();

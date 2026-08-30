(function () {
  const body = document.body;
  const toggle = document.querySelector("[data-nav-toggle]");
  const backdrop = document.querySelector("[data-nav-backdrop]");
  const filterInput = document.getElementById("nav-filter");
  const groups = Array.from(document.querySelectorAll(".nav-group"));

  function setNavOpen(open) {
    body.classList.toggle("nav-open", open);
    if (toggle) toggle.setAttribute("aria-expanded", open ? "true" : "false");
    if (backdrop) backdrop.hidden = !open;
  }

  if (toggle) {
    toggle.addEventListener("click", function () {
      setNavOpen(!body.classList.contains("nav-open"));
    });
  }

  if (backdrop) {
    backdrop.addEventListener("click", function () {
      setNavOpen(false);
    });
  }

  document.addEventListener("keydown", function (event) {
    if (event.key === "Escape") setNavOpen(false);
  });

  if (filterInput) {
    filterInput.addEventListener("input", function () {
      const query = filterInput.value.trim().toLowerCase();

      groups.forEach(function (group) {
        const links = Array.from(group.querySelectorAll(".nav-link"));
        let visible = 0;

        links.forEach(function (link) {
          const haystack = [
            link.dataset.title || "",
            link.dataset.type || "",
            group.dataset.topic || "",
          ]
            .join(" ")
            .toLowerCase();
          const match = !query || haystack.indexOf(query) !== -1;
          link.classList.toggle("is-hidden", !match);
          if (match) visible += 1;
        });

        group.classList.toggle("is-hidden", visible === 0);
      });
    });
  }

  const articleBody = document.querySelector("[data-article-body]");
  if (articleBody) {
    const firstBlock = articleBody.querySelector(":scope > blockquote");
    if (firstBlock && /^\s*Type\s*:/i.test(firstBlock.textContent || "")) {
      firstBlock.classList.add("type-source-tag");
    }
  }
})();

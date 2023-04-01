// Load the PDF file into the viewer
function loadPDF(url) {
    const pdfViewer = document.querySelector('.pdf-viewer');
    const pdfToolbar = document.querySelector('.pdf-toolbar');
    const closeBtn = document.querySelector('.close-pdf');
    const openExternalBtn = document.querySelector('.open-external');
    pdfjsLib.getDocument(url).promise.then((pdfDoc) => {
      const pdfViewer = document.querySelector('.pdf-viewer');
      const pdfPages = [];
      for (let i = 1; i <= pdfDoc.numPages; i++) {
        const pageDiv = document.createElement('div');
        pageDiv.classList.add('pdf-page');
        pdfPages.push(pageDiv);
        pdfViewer.appendChild(pageDiv);
        pdfDoc.getPage(i).then((page) => {
          const viewport = page.getViewport({ scale: 1 });
          const canvas = document.createElement('canvas');
          const ctx = canvas.getContext('2d');
          const renderContext = {
            canvasContext: ctx,
            viewport: viewport,
          };
          canvas.height = viewport.height;
          canvas.width = viewport.width;
          page.render(renderContext).promise.then(() => {
            pageDiv.appendChild(canvas);
          });
        });
      }
      pdfToolbar.style.display = 'flex';
      pdfViewer.style.display = 'block';
      closeBtn.addEventListener('click', () => {
        pdfViewer.innerHTML = '';
        pdfToolbar.style.display = 'none';
        pdfViewer.style.display = 'none';
      });
      openExternalBtn.addEventListener('click', () => {
        window.open(url, '_blank');
      });
    });
  }
  
  // Set up event listener for view PDF button
  const viewPdfBtns = document.querySelectorAll('.view-pdf');
  viewPdfBtns.forEach((btn) => {
    btn.addEventListener('click', (e) => {
      // Get PDF URL from button data attribute
      const url = e.target.dataset.url;
      loadPDF(url);
    });
  });
  
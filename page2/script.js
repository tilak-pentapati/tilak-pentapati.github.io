function loadPDF(url) {
    const pdfViewer = document.querySelector('.pdf-viewer');
    const pdfPages = document.querySelector('.pdf-pages');
    const pdfToolbar = document.querySelector('.pdf-toolbar');
    const closeBtn = document.querySelector('.close-pdf');
    const openExternalBtn = document.querySelector('.open-external');
    pdfjsLib.getDocument(url).promise.then((pdfDoc) => {
      for (let i = 1; i <= pdfDoc.numPages; i++) {
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
            pdfPages.appendChild(canvas);
          });
        });
      }
      pdfToolbar.style.display = 'flex';
      pdfViewer.style.display = 'block';
      closeBtn.addEventListener('click', () => {
        pdfPages.innerHTML = '';
        pdfToolbar.style.display = 'none';
        pdfViewer.style.display = 'none';
      });
      openExternalBtn.addEventListener('click', () => {
        window.open(url, '_blank');
      });
    });
  }
  
  window.addEventListener('load', () => {
    const pdfUrl = './cv.pdf';
    loadPDF(pdfUrl);
  });
  
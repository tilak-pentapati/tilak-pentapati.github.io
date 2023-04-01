function loadPDF(url) {
    const pdfPages = document.querySelector('#pdf-pages');
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
    });
  }
  
  window.addEventListener('load', () => {
    const pdfUrl = './CV_tilak.pdf';
    loadPDF(pdfUrl);
  });
  
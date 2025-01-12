const canvas = document.getElementById("canvas");
const canvasStream = canvas.captureStream(30);

const stetusElement = document.getElementById("movieStetus");

const mediaRecorder = new MediaRecorder(canvasStream);

canvas.addEventListener("startRender", () => {
    mediaRecorder.start();
    stetusElement.textContent = "：録画中";
});
canvas.addEventListener("finishRender", () => {
    mediaRecorder.stop();
    stetusElement.textContent = "：録画終了";
});

mediaRecorder.addEventListener("dataavailable", event => {
    const videoBlob = event.data;//new Blob([event.data], { type: event.data.type });
    const dataUrl = window.URL.createObjectURL(videoBlob);
    const anchor = document.getElementById("downloadLink");
    anchor.download = `disaster simulation movie`;
    anchor.href = dataUrl;
    stetusElement.textContent = "：準備完了";
});
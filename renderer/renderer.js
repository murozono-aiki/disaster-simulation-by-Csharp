import * as THREE from "three";
import { GLTFLoader } from "three/addons/loaders/GLTFLoader.js"

const width = 1280;
const height = 720;

// レンダラーを作成
const renderer = new THREE.WebGLRenderer({
    canvas: document.getElementById("canvas")
});
renderer.setSize(width, height);
renderer.setPixelRatio(window.devicePixelRatio);

// シーンを作成
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xFFFFFF);

// カメラを作成
const camera = new THREE.PerspectiveCamera(45, width / height, 1, 10000);
camera.position.set(300, 100, 1000);

// 平行光源
const light = new THREE.DirectionalLight(0xFFFFFF);
light.intensity = 2; // 光の強さを倍に
light.position.set(1, 1, 1); // ライトの方向
// 環境光
const ambientLight = new THREE.AmbientLight(0xbbbbbb);
scene.add(ambientLight);
// シーンに追加
scene.add(light);


let mouseX = 0; // マウス座標
let mouseY = 0; // マウス座標

/**
 * @param {MouseEvent} event
 */
function onMousemove_camera(event) {
    const currentMouseX = event.pageX;
    const currentMouseY = event.pageY;

    if (event.buttons % 2 == 1) {  // 左クリック時
        const cameraWorldDirection = new THREE.Vector3();
        camera.getWorldDirection(cameraWorldDirection);
        const cameraDirection = cameraWorldDirection.normalize();
        const cameraPositionVectorY = new THREE.Vector3(0, 1, 0);
        const cameraPositionVectorX = cameraDirection.clone().cross(cameraPositionVectorY.clone()).normalize();
        if (!event.ctrlKey) {
            const cameraMoveX = cameraPositionVectorX.clone().multiplyScalar(currentMouseX - mouseX);
            camera.position.x -= cameraMoveX.x;
            camera.position.y -= cameraMoveX.y;
            camera.position.z -= cameraMoveX.z;
            if (!event.shiftKey) {
                const cameraMoveY = cameraPositionVectorY.clone().multiplyScalar(currentMouseY - mouseY);
                camera.position.x += cameraMoveY.x;
                camera.position.y += cameraMoveY.y;
                camera.position.z += cameraMoveY.z;
            } else {
                const cameraMoveZ = cameraDirection.clone().multiplyScalar(currentMouseY - mouseY);
                camera.position.x += cameraMoveZ.x;
                camera.position.y += cameraMoveZ.y;
                camera.position.z += cameraMoveZ.z;
            }
        } else {
            if (Math.abs(currentMouseY - mouseY) > Math.abs(currentMouseX - mouseX)) {
                camera.rotateOnAxis(new THREE.Vector3(1, 0, 0), (currentMouseY - mouseY) / 80);
            } else {
                camera.rotateOnWorldAxis (new THREE.Vector3(0, 1, 0), (currentMouseX - mouseX) / 80);
            }
        }

        renderer.render(scene, camera);
    }

    mouseX = currentMouseX;
    mouseY = currentMouseY;
}

// マウス座標はマウスが動いた時のみ取得できる
document.addEventListener("mousemove", onMousemove_camera);


// 非同期処理で待機するのでasync function宣言とする
async function init() {
    // GLTF形式のモデルデータを読み込む
    const loader_building = new GLTFLoader();
    // GLTFファイルのパスを指定
    //const gltf = await loader.loadAsync('data/Juso-data.glb');
    const gltf_building = await loader_building.loadAsync('data/bldg_Building.glb');
    // 読み込み後に3D空間に追加
    const building = gltf_building.scene;
    building.position.set(0, 0, 0);
    scene.add(building);

    // GLTF形式のモデルデータを読み込む
    const loader_reliefFeature = new GLTFLoader();
    // GLTFファイルのパスを指定
    //const gltf = await loader.loadAsync('data/Juso-data.glb');
    const gltf_reliefFeature = await loader_reliefFeature.loadAsync('data/dem_ReliefFeature.glb');
    // 読み込み後に3D空間に追加
    const reliefFeature = gltf_reliefFeature.scene;
    reliefFeature.position.set(0, 0, 0);
    scene.add(reliefFeature);

    renderer.render( scene, camera );
}
let initPromise = init();


function startRecord() {
    const event = new Event("startRender");
    document.getElementById("canvas").dispatchEvent(event);
}
function finishRecord() {
    const event = new Event("finishRender");
    document.getElementById("canvas").dispatchEvent(event);
}



/** @type {{time:number, particlePositions:({x:number, y:number, z:number}|null)[]}[]} */
let resultData = [];
let currentIndex = 0;

async function loadResult() {
    const fileNameResponce = await fetch("./../source/output/fileNames.json");
    /** @type {string[]} */
    const fileNames = await fileNameResponce.json();
    for (let i = 0; i < fileNames.length; i++) {
        const resultResponce = await fetch("./../source/output/" + fileNames[i]);
        /** @type {{time:number, particlePositions:({x:number, y:number, z:number}|null)[]}[]} */
        const result = await resultResponce.json();
        resultData = resultData.concat(result);
    }
    renderParticles(resultData[0].particlePositions);
}
let loadPromise = loadResult();

function createParticleGeometry() {
    const geometry = new THREE.SphereGeometry(1, 3, 2);
    const material = new THREE.MeshStandardMaterial({color: 0x0000FF});
    const sphere = new THREE.Mesh(geometry, material);
    scene.add(sphere);
    return sphere;
}
/** @type {THREE.Mesh<THREE.SphereGeometry, THREE.MeshStandardMaterial, THREE.Object3DEventMap>[]} */
const particleGeometries = [];

/** @param {({x:number, y:number, z:number}|null)[]} particlePositions */
function renderParticles(particlePositions) {
    for (let i = 0; i < Math.max(particlePositions.length, particleGeometries.length); i += 1) {
        if (particlePositions[i]) {
            if (!particleGeometries[i]) {
                particleGeometries[i] = createParticleGeometry();
            }
            particleGeometries[i].position.x = particlePositions[i].x;
            particleGeometries[i].position.y = particlePositions[i].y;
            particleGeometries[i].position.z = particlePositions[i].z;
        } else {
            if (particleGeometries[i]) {
                particleGeometries[i].geometry.dispose();
                particleGeometries[i] = undefined;
            }
        }
    }
    renderer.render(scene, camera);
}
let start;
function renderStep(timestamp) {
    if (start === undefined) {
        startRecord();
        start = timestamp;
    }
    const elapsed = (timestamp - start) / 1000;  // 秒単位

    let currentResultData;
    while (currentIndex < resultData.length && (!resultData[currentIndex] || resultData[currentIndex].time <= elapsed)) {
        currentIndex++;
        if (resultData[currentIndex]) currentResultData = resultData[currentIndex];
    }
    if (currentResultData) {
        document.getElementById("time").textContent = Math.round(currentResultData.time) + "秒";
        renderParticles(currentResultData.particlePositions);
    }

    if (currentIndex < resultData.length) {
        window.requestAnimationFrame(renderStep);
    } else {
        finishRecord();
    }
}

Promise.all([initPromise, loadPromise]).then(() => {
    /** @type {HTMLButtonElement} */
    const startButton = document.getElementById("start");
    startButton.disabled = false;
    startButton.addEventListener("click", event => {
        currentIndex = 0;
        window.requestAnimationFrame(renderStep);
    });
});
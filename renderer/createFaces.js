import * as THREE from "three";
import { GLTFLoader } from "three/addons/loaders/GLTFLoader.js";
import { BlobReader, BlobWriter, TextReader, ZipWriter } from "./zip.js-2.7.54/index.js";

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

    const normalVectors_building = createNormalVectors(building.children[0].geometry);
    const normalVectors_reliefFeature = createNormalVectors(reliefFeature.children[0].geometry);

    const normalVectors = [].concat(normalVectors_building, normalVectors_reliefFeature);
    const normalVectors_sortX = normalVectors.toSorted((a, b) => a.centerOfGravity.x - b.centerOfGravity.x);
    const normalVectors_sortY = normalVectors.toSorted((a, b) => a.centerOfGravity.y - b.centerOfGravity.y);
    const normalVectors_sortZ = normalVectors.toSorted((a, b) => a.centerOfGravity.z - b.centerOfGravity.z);

    const particle = createParticleGeometry();
    particle.position.x = -700;
    particle.position.y = 45;
    particle.position.z = 2303;

    const zipFileWriter = new BlobWriter('application/json');
    const textReaders = [];
    while (normalVectors.length > 0) {
        const array = normalVectors.splice(0, 70000);
        const json = JSON.stringify(array);
        const textReader = new TextReader(json);
        textReaders.push(textReader);
    }
    const zipWriter = new ZipWriter(zipFileWriter);
    for (let i = 0; i < textReaders.length; i++) {
        await zipWriter.add("face-data-" + (i + 1) + ".json", textReaders[i]);
    }
    await zipWriter.close();
    const zipFileBlob = await zipFileWriter.getData();
    const url = URL.createObjectURL(zipFileBlob);
    const anchor = document.getElementById("downloadLink");
    anchor.download = `face-data.zip`;
    anchor.href = url;
    anchor.textContent = "面データをダウンロード";
}
init();

function createNormalVectors(geometry) {
    const result = [];

    const positionAttributes = geometry.attributes.position;
    const indexArray = geometry.index.array;

    for (let i = 0; i < indexArray.length; i += 3) {
        const point1 = new THREE.Vector3(positionAttributes.getX(indexArray[i]), positionAttributes.getY(indexArray[i]), positionAttributes.getZ(indexArray[i]));
        const point2 = new THREE.Vector3(positionAttributes.getX(indexArray[i + 1]), positionAttributes.getY(indexArray[i + 1]), positionAttributes.getZ(indexArray[i + 1]));
        const point3 = new THREE.Vector3(positionAttributes.getX(indexArray[i + 2]), positionAttributes.getY(indexArray[i + 2]), positionAttributes.getZ(indexArray[i + 2]));

        const centerOfGravity = point1.clone().add(point2.clone()).add(point3.clone()).divideScalar(3);

        const normalVector = point2.clone().sub(point1.clone()).cross(point3.clone().sub(point1.clone())).normalize();
        if (normalVector.x == 0 && normalVector.y == 0 && normalVector.z == 0) {
            continue;
        }

        const insideJudge = [];
        const dampingCoefficient = -5;  // ダンパ係数
        const height = 5;  // 最大食い込み量（食い込んでいる面を特定するために用いる）←他の方法が望ましいが、簡略化のため
        insideJudge.push({
            point: point1,
            normalVector: multiplyScalarVector3(normalVector, -1)
        });
        insideJudge.push({
            point: subVector3(point1, multiplyScalarVector3(normalVector, height)),
            normalVector: normalVector
        });
        /*insideJudge.push({
            point: addVector3(point1, multiplyScalarVector3(normalVector, -dampingCoefficient)),
            normalVector: multiplyScalarVector3(normalVector, -1)
        });
        insideJudge.push({
            point: subVector3(point1, multiplyScalarVector3(normalVector, -dampingCoefficient)),
            normalVector: normalVector
        });*/
        insideJudge.push({
            point: point1,
            normalVector: normalVector.clone().cross(point2.clone().sub(point1.clone())).normalize()
        });
        insideJudge.push({
            point: point2,
            normalVector: normalVector.clone().cross(point3.clone().sub(point2.clone())).normalize()
        });
        insideJudge.push({
            point: point3,
            normalVector: normalVector.clone().cross(point1.clone().sub(point3.clone())).normalize()
        });

        result.push({
            index: i / 3,
            triangle: [point1, point2, point3],
            centerOfGravity: centerOfGravity,
            normalVector: normalVector,
            insideJudge: insideJudge
        });
    }
    return result;
}

function createParticleGeometry() {
    const geometry = new THREE.SphereGeometry(3, 3, 2);
    const material = new THREE.MeshStandardMaterial({color: 0x0000FF});
    const sphere = new THREE.Mesh(geometry, material);
    scene.add(sphere);
    return sphere;
}

/**
 * ベクトルについて a + b を計算する関数
 * @param {Vector3} a - 足すベクトル
 * @param {Vector3} b - 足されるベクトル
 * @returns {Vector3} - ベクトルの和
 */
function addVector3(a, b) {
    return createVector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

/**
 * ベクトルについて a - b を計算する関数
 * @param {Vector3} a - 引かれるベクトル
 * @param {Vector3} b - 引くベクトル
 * @returns {Vector3} - ベクトルの差
 */
function subVector3(a, b) {
    return createVector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

/**
 * ベクトルについて aのn倍 を計算する関数
 * @param {Vector3} a - ベクトル
 * @param {number} n - 実数（スカラー）
 * @returns {Vector3} - 実数倍したベクトル
 */
function multiplyScalarVector3(a, n) {
    return createVector3(n * a.x, n * a.y, n * a.z);
}
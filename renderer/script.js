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

    const normalVectors_building = createNormalVectors(building.children[0].geometry);
    const normalVectors_reliefFeature = createNormalVectors(reliefFeature.children[0].geometry);

    const normalVectors = [].concat(normalVectors_building, normalVectors_reliefFeature);
    const normalVectors_sortX = normalVectors.toSorted((a, b) => a.centerOfGravity.x - b.centerOfGravity.x);
    const normalVectors_sortY = normalVectors.toSorted((a, b) => a.centerOfGravity.y - b.centerOfGravity.y);
    const normalVectors_sortZ = normalVectors.toSorted((a, b) => a.centerOfGravity.z - b.centerOfGravity.z);

    /*const material = new THREE.MeshStandardMaterial({color: 0x0000FF});
    for (let i = 655000; i < 700000; i++) {
        const normalVectorObj = normalVectors[i];
        const vector = normalVectorObj.normalVector.multiplyScalar(5);
        for (let j = 0; j < normalVectorObj.triangle.length; j++) {
            const points = [];
            points.push(normalVectorObj.triangle[j].clone());
            points.push(normalVectorObj.triangle[j].clone().add(vector.clone()));

            const geometry = new THREE.BufferGeometry().setFromPoints( points );
            const line = new THREE.Line( geometry, material );
            scene.add( line );
        }
    }*/

    renderer.render( scene, camera );

    simulateStart({
        normalVectors: {
            data: normalVectors,
            sort: {
                x: normalVectors_sortX,
                y: normalVectors_sortY,
                z: normalVectors_sortZ
            }
        }
    });
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


function startRecord() {
    const event = new Event("startRender");
    document.getElementById("canvas").dispatchEvent(event);
}
function finishRecord() {
    const event = new Event("finishRender");
    document.getElementById("canvas").dispatchEvent(event);
}



function createParticleGeometry() {
    const geometry = new THREE.SphereGeometry(1, 3, 2);
    const material = new THREE.MeshStandardMaterial({color: 0x0000FF});
    const sphere = new THREE.Mesh(geometry, material);
    scene.add(sphere);
    return sphere;
}
const particleGeometries = [];

function simulateStart(data) {
    const simulateWorker = new Worker("simulator.js");
    simulateWorker.postMessage(data);
    let progressReported = false;
    simulateWorker.addEventListener("message", event => {
        const object = event.data;
        switch (object.type) {
            case "result":
                const particles = object.content;
                for (let i = 0; i < Math.max(particles.length, particleGeometries.length); i += 1) {
                    if (particles[i] && !particles[i].is_wall) {
                        if (!particleGeometries[i]) {
                            particleGeometries[i] = createParticleGeometry();
                        }
                        particleGeometries[i].position.x = particles[i].position.x;
                        particleGeometries[i].position.y = particles[i].position.y;
                        particleGeometries[i].position.z = particles[i].position.z;
                    } else {
                        if (particleGeometries[i]) {
                            particleGeometries[i].dispose();
                            particleGeometries[i] = undefined;
                        }
                    }
                }
                renderer.render(scene, camera);
                break;
            case "progress":
                if (!progressReported) {
                    progressReported = true;
                    requestAnimationFrame(() => {
                        document.getElementById("progress").textContent = object.content;
                    });
                    setTimeout(() => {
                        progressReported = false;
                    }, 5000);
                }
                break;
            case "sphere":
                const particle = createParticleGeometry();
                particle.position.x = object.point.x;
                particle.position.y = object.point.y;
                particle.position.z = object.point.z;
        }
    });
}




/* シミュレーション部分 */

/**
 * 加速度を計算するときに使う項
 * @type {{pressureTerm:THREE.Vector3, viscosityTerm:THREE.Vector3}[]}
 */
const terms = [];

const h = 0.012; //影響半径
const particleMass = 0.0002; //粒子の質量

const g = new THREE.Vector3(0, -9.8, 0);  // 重力加速度
const pressureStiffness = 200; //圧力係数
const restDensity = 1000; //静止密度
const viscosity = 1;  // 粘性係数

const densityCoef = particleMass * 315 / (64 * Math.PI * Math.pow(h,9)); //密度計算で使うヤツ

const pressureCoef = particleMass * 45 / (Math.PI * Math.pow(h,6)); //圧力項計算で使うヤツ
const viscosityCoef = viscosity * particleMass * 45 / (Math.PI * Math.pow(h,6)); //粘性項計算で使うヤツ

//壁は立方体
const wallWidth = 10; 
const thickness = 0.1/* 0.5 */; //壁となる粒子の厚み（個数）
const interval = h / thickness;


const shaft = ["x","y","z"];

//壁となる、速度がゼロで固定の粒子群
function makeWall(particles){
    let newPosition = new THREE.Vector3(0,0,0);
    for(let x = 0; x < 11 ; x = x + 10){
        newPosition[shaft[2]] = x;
        for(let i=0 ; i< wallWidth / interval; i++){//高さ
            for(let j = 0; j < wallWidth / interval; j++){//横幅
                newPosition[shaft[1]] = (interval / 2) + i * interval;
                newPosition[shaft[0]] = (interval / 2) + j * interval;
                let newParticle = {
                    position: newPosition.clone(),
                    velocity: new THREE.Vector3(0, 0, 0),
                    force: new THREE.Vector3(0, 0, 0),
                    density: 0,
                    pressure: 0,
                    is_wall: true
                };
                particles.push(newParticle);
            }

        }  
    }
    for(let x = 0; x < 11 ; x = x + 10){
        newPosition[shaft[0]] = x;
        for(let i=0 ; i< wallWidth / interval; i++){//高さ
            for(let j = 0; j < wallWidth / interval; j++){//横幅
                newPosition[shaft[1]] = (interval / 2) + i * interval;
                newPosition[shaft[2]] = (interval / 2) + j * interval;
                let newParticle = {
                    position: newPosition.clone(),
                    velocity: new THREE.Vector3(0, 0, 0),
                    force: new THREE.Vector3(0, 0, 0),
                    density: 0,
                    pressure: 0,
                    is_wall: true
                };
                particles.push(newParticle);
            }

        }  
    }
    newPosition[shaft[1]] = 0;
        for(let i=0 ; i< wallWidth / interval; i++){//x方向
            for(let j = 0; j < wallWidth / interval; j++){//z方向
                newPosition[shaft[0]] = (interval / 2) + i * interval;
                newPosition[shaft[2]] = (interval / 2) + j * interval;
                let newParticle = {
                    position: newPosition.clone(),
                    velocity: new THREE.Vector3(0, 0, 0),
                    force: new THREE.Vector3(0, 0, 0),
                    density: 0,
                    pressure: 0,
                    is_wall: true
                };
                particles.push(newParticle);
            }

        }
    console.log(particles)
}
 
/**
 * 粒子の密度計算
 * @param {{position:THREE.Vector3, velocity:THREE.Vector3, force:THREE.Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcDensity(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算
        let nowParticle = particles[i]; //今回計算する粒子
        let sum = 0; //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = nearParticle.position.clone().sub(nowParticle.position.clone()); //粒子距離
            const r2 = diff.clone().dot(diff.clone()); //粒子距離の２乗
 
            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const c = h2 - r2;
                sum += Math.pow(c,3); //(h2-r2)の３乗
            }
        }
 
        nowParticle.density = sum * densityCoef; //密度が求まった
    }
}

/**
 * 粒子の圧力計算
 * @param {{position:THREE.Vector3, velocity:THREE.Vector3, force:THREE.Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcPressure(particles) {
    for(let i = 0; i < particles.length; i++) { //一つづつ粒子の圧力を計算
        particles[i].pressure = pressureStiffness * ( particles[i].density - restDensity );
    }
}

/**
 * 粒子の圧力項計算
 * @param {{position:THREE.Vector3, velocity:THREE.Vector3, force:THREE.Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcPressureTerm(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算
        let nowParticle = particles[i]; //今回計算する粒子
        let sum = new THREE.Vector3(0, 0, 0); //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = nearParticle.position.clone().sub(nowParticle.position.clone()); //粒子距離
            const r2 = diff.clone().dot(diff.clone()); //粒子距離の２乗
 
            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const r = Math.sqrt(r2); //粒子距離
                const c = h - r;
                const n = ((nearParticle.pressure /*-*/+ nowParticle.pressure) / (2 * nearParticle.density)) * Math.pow(c,2) / r;
                sum = sum.add(diff.multiplyScalar(n));
            }
        }
 
        if (!terms[i]) terms[i] = {};
        terms[i].pressureTerm = sum.multiplyScalar((-1/*/nowParticle.pressure*/) * pressureCoef);  // 圧力項が求まった
    }
}

/**
 * 粒子の粘性項計算
 * @param {{position:THREE.Vector3, velocity:THREE.Vector3, force:THREE.Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function calcViscosityTerm(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算
        let nowParticle = particles[i]; //今回計算する粒子
        let sum = new THREE.Vector3(0, 0, 0); //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = nearParticle.position.clone().sub(nowParticle.position.clone()); //粒子距離
            const r2 = diff.clone().dot(diff.clone()); //粒子距離の２乗
 
            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const r = Math.sqrt(r2); //粒子距離
                const c = h - r;
                const n = c / nearParticle.density;
                sum = sum.add(nearParticle.velocity.clone().sub(nowParticle.velocity.clone()).multiplyScalar(n));
            }
        }
 
        if (!terms[i]) terms[i] = {};
        terms[i].viscosityTerm = sum.multiplyScalar(viscosityCoef);  // 粘性項が求まった
    }
}


/**
 * 粒子のリスト
 * @type {{
 *      position:THREE.Vector3,
 *      velocity:THREE.Vector3,
 *      force:THREE.Vector3,
 *      density:number,
 *      pressure:number,
 *      is_wall:boolean
 * }[]}
 */
const _particles = [];

async function start() {
    makeWall(_particles);
}
//start();

const defaultDeltaTime = 0.1;

let previousTimeStamp;
async function tick(timestamp) {
    console.log("tick");
    calcDensity(_particles);
    calcPressure(_particles);
    calcPressureTerm(_particles);
    calcViscosityTerm(_particles);

    if (previousTimeStamp === undefined) previousTimeStamp = timestamp;
    const deltaTime = timestamp - previousTimeStamp;
    for (let i = 0; i < _particles.length; i++) {
        const nowParticle = _particles[i];
        const a = terms[i].pressureTerm.add(terms[i].viscosityTerm).add(g);
        const v = nowParticle.velocity.clone().add(a.multiplyScalar(deltaTime));
        nowParticle.velocity = v;
    }
    

    // レンダリング
    renderer.render(scene, camera);

    previousTimeStamp = timestamp;
    if (previousTimeStamp < 0.5) {
        tick(previousTimeStamp + defaultDeltaTime);
    }
}

// 初回実行
//tick(0);
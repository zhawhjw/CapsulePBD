import * as THREE from "three";
import * as PHY from "simplePhysics";
import {OrbitControls} from "three/addons/controls/OrbitControls.js";

import Stats from "three/addons/libs/stats.module.js";

let renderer, scene, camera;
let world = {
  x: 80,
  z: 80,
};
let agentData = [];
let pickableObjects = [];
let selected = null;
let mouse = new THREE.Vector2();
const raycaster = new THREE.Raycaster();
let grid, ring;

let spotLights = {};
let topTextures = {};
let topTexture;
const RADIUS = 1;
const blueAgentMaterial = new THREE.MeshLambertMaterial({
  color: 0x0000ff,
});
const redAgentMaterial = new THREE.MeshLambertMaterial({
  color: 0xff0000,
});
const greenAgentMaterial = new THREE.MeshLambertMaterial({
  color: 0x00ff00,
});

const stats = new Stats();
document.body.appendChild(stats.dom);

init();
render();

function init() {
  // renderer
  renderer = new THREE.WebGLRenderer();
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap; //
  renderer.setSize(window.innerWidth, window.innerHeight);
  document.body.appendChild(renderer.domElement);

  // scene
  scene = new THREE.Scene();
  // camera
  camera = new THREE.PerspectiveCamera(
    45,
    window.innerWidth / window.innerHeight,
    1,
    1000
  );

  camera.position.set(-67.26, 54.07, -3.77);
  camera.rotation.order = "YXZ";
  camera.rotation.y = -1.6267;
  camera.rotation.x = -0.46;

  // controls
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.addEventListener("change", render);
  controls.enableZoom = false;
  controls.enablePan = false;
  controls.maxPolarAngle = Math.PI / 2;

  // light
  const light = new THREE.PointLight(0xffffff, 0.9, 0, 100000);
  light.position.set(0, 50, 120);
  light.castShadow = true;
  light.shadow.mapSize.width = 512; // default
  light.shadow.mapSize.height = 512; // default
  light.shadow.camera.near = 0.5; // default
  light.shadow.camera.far = 5000; // default

  const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
  directionalLight.castShadow = true;
  directionalLight.position.set(-5, 20, 4);
  directionalLight.target.position.set(9, 0, -9);
  directionalLight.shadow.camera.left *= 9;
  directionalLight.shadow.camera.right *= 9;
  directionalLight.shadow.camera.top *= 9;
  directionalLight.shadow.camera.bottom *= 9;

  scene.add(directionalLight);

  // axes
  scene.add(new THREE.AxesHelper(40));
  const loader = new THREE.TextureLoader();
  const texture = loader.load("resources/OIP.jpg");
  texture.wrapS = THREE.RepeatWrapping;
  texture.wrapT = THREE.RepeatWrapping;
  texture.magFilter = THREE.NearestFilter;
  const repeats = 40 / 32;
  texture.repeat.set(repeats, repeats);

  topTexture = loader.load("resources/triangle2.png");
  //topTexture.wrapS = THREE.RepeatWrapping;
  //topTexture.wrapT = THREE.RepeatWrapping;
  topTexture.magFilter = THREE.NearestFilter;
  topTexture.repeat.set(3, 3);
  //topTexture.rotation = -Math.PI / 2;
  // grid
  const geometry = new THREE.PlaneGeometry(world.x, world.z, 10, 10);
  const material = new THREE.MeshPhongMaterial({
    map: texture,
    //side: THREE.DoubleSide,
  });
  grid = new THREE.Mesh(geometry, material);
  grid.castShadow = true; //default is false
  grid.receiveShadow = true; //default
  grid.rotation.order = "YXZ";
  grid.rotation.y = -Math.PI / 2;
  grid.rotation.x = -Math.PI / 2;
  scene.add(grid);

  // experiment border
  // const boxGeometry1 = new THREE.BoxGeometry(50, 5, 1);
  // const boxMaterial1 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
  // const left = new THREE.Mesh(boxGeometry1, boxMaterial1);
  // left.position.set(0, 2.5, -25);
  //
  // const boxGeometry2 = new THREE.BoxGeometry(1, 5, 50);
  // const boxMaterial2 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
  // const bottom = new THREE.Mesh(boxGeometry2, boxMaterial2);
  // bottom.position.set(-25, 2.5, 0);
  //
  // const boxGeometry3 = new THREE.BoxGeometry(1, 5, 50);
  // const boxMaterial3 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
  // const top = new THREE.Mesh(boxGeometry3, boxMaterial3);
  // top.position.set(25, 2.5, 0);
  //
  // const boxGeometry4 = new THREE.BoxGeometry(50, 5, 1);
  // const boxMaterial4 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
  // const right = new THREE.Mesh(boxGeometry4, boxMaterial4);
  // right.position.set(0, 2.5, 25);
  //
  // scene.add(left);
  // scene.add(bottom);
  // scene.add(top);
  // scene.add(right);

  const ringGeometry = new THREE.RingGeometry(1, 3, 12);
  const ringMaterial = new THREE.MeshBasicMaterial({
    color: 0xffff00,
    side: THREE.DoubleSide,
  });
  ring = new THREE.Mesh(ringGeometry, ringMaterial);
  scene.add(ring);
  ring.rotation.x = -Math.PI / 2;
  ring.position.y += 0.01;


  function sampleCirclePoints(radius, sampleCount, centerX, centerY) {
    let points = [];
    for (let i = 0; i < sampleCount; i++) {
      // Angle in radians
      let angle = 2 * Math.PI * i / sampleCount;

      // Calculating x and y coordinates
      let x = centerX + radius * Math.cos(angle);
      let y = centerY + radius * Math.sin(angle);

      points.push({ x: x, y: y });
    }
    return points;
  }


  function sampleCirclePointsWithDistance(sampleCount, distanceBetweenPoints, centerX, centerY) {
    // Calculating the total circumference that would fit the points with the given distance
    let totalCircumference = distanceBetweenPoints * sampleCount;

    // Updating the radius based on the new circumference
    let radius = totalCircumference / (2 * Math.PI);

    let points = [];
    for (let i = 0; i < sampleCount; i++) {
      // Angle in radians
      let angle = 2 * Math.PI * i / sampleCount;

      // Calculating x and y coordinates
      let x = centerX + radius * Math.cos(angle);
      let y = centerY + radius * Math.sin(angle);

      points.push({ x: x, y: y });
    }
    return points;
  }



  function testScenario(){
    addColumnAgentGroup(
        agentData,
        1,
        RADIUS * 4,
        {
          x: 20,
          z: 0,
        },
        {
          x: -20,
          z: -5,
        },
        0.8,
        "X"
    );

    addColumnAgentGroup(
        agentData,
        1,
        RADIUS * 4,
        {
          x: -20,
          z: 0,
        },
        {
          x: 20,
          z: 0,
        },
        0.8,
        "X"
    );
  }

  function testHallwayScenario(){

    for(let i=0;i<1;i++){
      addColumnAgentGroup(
          agentData,
          1,
          RADIUS * 1.5,
          {
            x: 20 ,
            z: 0 + i * 6,
          },
          {
            x: -20,
            z: 0 + i * 6,
          },
          0.8,
          "X"
      );

      addColumnAgentGroup(
          agentData,
          1,
          RADIUS * 4,
          {
            x: -20,
            z: 0 + i * 6 - 2
          },
          {
            x: 20,
            z: 0 + i * 6  - 2
          },
          0.8,
          "X"
      );


    }


  }

  function testCrossScenario(){
    addColumnAgentGroup(
        agentData,
        1,
        RADIUS * 4,
        {
          x: 20,
          z: 0,
        },
        {
          x: -20,
          z: -5,
        },
        0.8,
        "X"
    );

    addColumnAgentGroup(
        agentData,
        1,
        RADIUS * 4,
        {
          x: 0,
          z: 20,
        },
        {
          x: 0,
          z: -20,
        },
        0.8,
        "X"
    );
  }
  function testCrossWithDiagnoScenario(){
    addColumnAgentGroup(
        agentData,
        1,
        RADIUS * 4,
        {
          x: 20,
          z: 0,
        },
        {
          x: -20,
          z: -5,
        },
        0.8,
        "X"
    );

    addColumnAgentGroup(
        agentData,
        1,
        RADIUS * 4,
        {
          x: -10,
          z: 10,
        },
        {
          x: 20,
          z: -20,
        },
        0.8,
        "X"
    );
  }

  function defaultScenario(){
    addColumnAgentGroup(
        agentData,
        50,
        RADIUS * 4,
        {
          x: 0,
          z: 0,
        },
        {
          x: -10,
          z: 10,
        },
        0.8,
        "X"
    );

  }

  function circleScenario(){

    let points = sampleCirclePoints(20, 20, 0, 0);
    // let points = sampleCirclePointsWithDistance(42, 2 * 2* RADIUS + 2, 0, 0);
    console.log(points);
    points.forEach(function (point){
      addColumnAgentGroup(
          agentData,
          1,
          0,
          {
            x: point.x,
            z: point.y,
          },
          {
            x: -point.x ,
            z: -point.y,
          },
          0.8,
          "X"
      );
    });

  }

  function addColumnAgentGroup(
    agentData,
    numAgents,
    spacing,
    startPos,
    goalPos,
    velocityMagnitude,
    direction
  ) {
    let i = 0;
    let initalIdx = agentData.length;
    let dx = 0,
      dz = 0,
      vx = 0,
      vz = 0;
    let distanceToGoal = PHY.distance(
      startPos.x,
      startPos.z,
      goalPos.x,
      goalPos.z
    );
    vx = (velocityMagnitude * (goalPos.x - startPos.x)) / distanceToGoal;
    vz = (velocityMagnitude * (goalPos.z - startPos.z)) / distanceToGoal;

    if (direction == "X") {
      dx = spacing;
    } else if (direction == "Z") {
      dz = spacing;
    }
    while (i < numAgents) {
      agentData.push({
        index: i + initalIdx,
        x: startPos.x + dx * i,
        y: 2.0,
        z: startPos.z + dz * i,
        goal_x: goalPos.x + dx * i,
        goal_y: 0.0,
        goal_z: goalPos.z + dz * i,
        vx: vx,
        vy: 0.0,
        vz: vz,
        v_pref: Math.sqrt(vx * vx + vz * vz),
        radius: RADIUS,
        invmass: 0.5,
        colliding: false,
        group_id: 1,

        best: null
      });
      i += 1;
    }
  }

  let i = 0;
  let deltaSpacing = 3;
  let startX, startY, goalX, goalY;
  startX = -25;
  goalX = -25;
  startY = -20;
  goalY = 20;
  world.distanceConstraints = [];

  // defaultScenario();
  // testScenario()
  testCrossScenario();
  // testCrossWithDiagnoScenario();
  // testHallwayScenario();
  // testCrossScenario();
  // circleScenario();



  let agentGeom, agentMaterial, agent;
  let spotLight, spotLightTarget;
  let agentPointGeom, agentPointMaterial, agentPoint;

  agentData.forEach(function (item) {
    //agentGeom = new THREE.CylinderGeometry(item.radius, 1, 4, 16);
    agentGeom = new THREE.CapsuleGeometry(item.radius, 2 * item.radius, 4, 8);
    agentMaterial = new THREE.MeshLambertMaterial({
      color: 0x00ff00,
    });
    agent = new THREE.Mesh(agentGeom, agentMaterial);
    agent.castShadow = true;
    agent.receiveShadow = true;
    agent.userData = {
      index: item.index,
    };
    agent.rotateX(Math.PI / 2);
    // agent.rotateZ(Math.PI / 2);
    scene.add(agent);



    agentPointGeom = new THREE.CapsuleGeometry(item.radius, 2 * item.radius, 4, 8);
    agentPointMaterial = new THREE.MeshLambertMaterial({
      color: 0x0000ff,
    });
    agentPoint = new THREE.Mesh(agentPointGeom, agentPointMaterial);
    agentPoint.castShadow = true;
    agentPoint.receiveShadow = true;
    scene.add(agentPoint);




    // -----------------
    //adding spotlight code
    spotLight = new THREE.SpotLight(0xffffff);
    spotLight.position.set(item.x, item.y + 6, item.z);
    spotLight.shadow.mapSize.width = 1024;
    spotLight.shadow.mapSize.height = 1024;
    spotLight.shadow.camera.near = 500;
    spotLight.shadow.camera.far = 4000;
    spotLight.shadow.camera.fov = 30;
    spotLight.intensity = 0.4;
    spotLight.angle = Math.PI / 8;
    spotLightTarget = new THREE.Object3D();
    scene.add(spotLightTarget);
    spotLight.target = spotLightTarget;
    scene.add(spotLight);
    spotLights[item.index] = spotLight;
    // ----------------
    item.agent = agent;
    item.agentPoint = agentPoint;
    pickableObjects.push(agent);
    // pickableObjects.push(agent);

  });
  window.addEventListener("resize", onWindowResize);
  window.addEventListener("mousedown", mouseDown, false);
  window.addEventListener("mousemove", mouseMove, false);
}

function onWindowResize() {
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
}

function mouseMove(event) {
  event.preventDefault();
  if (selected != null) {
    mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
    raycaster.setFromCamera(mouse, camera);
    var intersects = raycaster.intersectObject(grid, false);
    for (let i = 0; i < intersects.length; i++) {

      agentData.forEach(function (member) {
        if (selected != null && member.index === selected) {
          member.goal_x = intersects[i].point.x;
          member.goal_z = intersects[i].point.z;
          // ring.position.x = intersects[i].point.x;
          // ring.position.z = intersects[i].point.z;
        }
      });
      break;
    }
  }
}

function mouseDown(event) {
  mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
  mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
  raycaster.setFromCamera(mouse, camera);
  selected = null;
  var intersects = raycaster.intersectObjects(pickableObjects, false);
  for (var i = 0; i < intersects.length; i++) {
    /* TODO finish this part as
     */
    selected = intersects[i].object.userData.index;
    break;
  }
}

function render() {
  renderer.render(scene, camera);
}

function animate() {
  // console.log(mouse.x, mouse.y);
  requestAnimationFrame(animate);
  PHY.step(RADIUS, agentData, world, scene);



  agentData.forEach(function (member) {
    // prevent agents from leaving the walls
    // if (member.x > 23) {
    //   member.x = 23;
    // }
    // if (member.x < -23) {
    //   member.x = -23;
    // }
    // if (member.z > 23) {
    //   member.z = 23;
    // }
    // if (member.z < -23) {
    //   member.z = -23;
    // }

    member.agent.position.x = member.x;
    member.agent.position.y = member.y;
    member.agent.position.z = member.z;

    if(member.best !== null){
      member.agentPoint.position.x = member.best.x;
      member.agentPoint.position.y = member.y;
      member.agentPoint.position.z = member.best.z;
    }


    const dx = member.goal_x - member.x;
    const dz = member.goal_z - member.z;
    member.agent.rotation.z = Math.atan2(dz, dx);

    member.agent.material = redAgentMaterial;
    if (member.colliding) {
      member.agent.material = greenAgentMaterial;
    }
    member.colliding = false;

    if (selected != null && member.index === selected) {
      member.agent.material = blueAgentMaterial;
    }
    /* TODO finish this part as 
        

         */
    spotLights[member.index].position.set(
      member.x - member.vx,
      member.y - member.vy,
      member.z - member.vz
    );
    spotLights[member.index].target.position.x = member.x;
    spotLights[member.index].target.position.y = member.y;
    spotLights[member.index].target.position.z = member.z;
  });
  renderer.render(scene, camera);
  stats.update();
}

animate();

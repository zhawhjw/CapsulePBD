import * as THREE from "three";

export function distance(x1, y1, x2, y2) {
  return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

export function step(RADIUS, sceneEntities, world, scene, customParams = {}) {
  const AGENTSIZE = RADIUS * 2;
  const epsilon = 0.0001;
  const timestep = 0.03;
  const ITERNUM = 1; // 3
  const agentLength = RADIUS;
  const MAX_DELTA=AGENTSIZE;

  // collision functions

  function rotateLineSegment(x1, y1, x2, y2, r) {
    // Calculate the center of the line segment
    const centerX = (x1 + x2) / 2;
    const centerY = (y1 + y2) / 2;

    // Translate the line segment so that its center is at the origin
    const x1p = x1 - centerX;
    const y1p = y1 - centerY;
    const x2p = x2 - centerX;
    const y2p = y2 - centerY;

    // Rotate the line segment about the origin
    const cosR = Math.cos(r);
    const sinR = Math.sin(r);
    const x1r = x1p * cosR - y1p * sinR;
    const y1r = x1p * sinR + y1p * cosR;
    const x2r = x2p * cosR - y2p * sinR;
    const y2r = x2p * sinR + y2p * cosR;

    // Translate the line segment back to its original position
    const newX1 = x1r + centerX;
    const newY1 = y1r + centerY;
    const newX2 = x2r + centerX;
    const newY2 = y2r + centerY;

    // Return the new endpoints of the line segment
    return [newX1, newY1, newX2, newY2];
  }

  function ClosestPointOnLineSegment(A, B, Point) {
    const AB = B.clone().sub(A);
    const t = AB.clone().dot(Point.clone().sub(A)) / AB.clone().dot(AB);
    return A.clone().add(
      AB.clone().multiplyScalar(Math.min(Math.max(t, 0), 1))
    );
  }

  function PointOnLineSegment(A, B, Point) {
    const AB = B.clone().sub(A);
    const t = AB.clone().dot(Point.clone().sub(A)) / AB.clone().dot(AB);
    return A.clone().add(
        AB.clone().multiplyScalar(t)
    );
  }

  function is_colliding_torso(x11, y11, x12, y12, x21, y21, x22, y22) {
    // console.log(segments_distance(x11, y11, x12, y12, x21, y21, x22, y22));
    return (
      segments_distance(x11, y11, x12, y12, x21, y21, x22, y22) < 2 * RADIUS
    );
  }

  function segments_distance(x11, y11, x12, y12, x21, y21, x22, y22) {
    // distance between two segments in the plane:
    // one segment is (x11, y11) to (x12, y12)
    // the other is (x21, y21) to (x22, y22)

    if (segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22)) return 0;

    // try each of the 4 vertices w/the other segment
    let distances = [
      point_segment_distance(x11, y11, x21, y21, x22, y22),
      point_segment_distance(x12, y12, x21, y21, x22, y22),
      point_segment_distance(x21, y21, x11, y11, x12, y12),
      point_segment_distance(x22, y22, x11, y11, x12, y12),
    ];
    return Math.min(...distances);
  }

  function segments_intersect(x11, y11, x12, y12, x21, y21, x22, y22) {
    // whether two segments in the plane intersect:
    // one segment is (x11, y11) to (x12, y12)
    // the other is (x21, y21) to (x22, y22)

    let dx1 = x12 - x11;
    let dy1 = y12 - y11;
    let dx2 = x22 - x21;
    let dy2 = y22 - y21;

    let delta = dx2 * dy1 - dy2 * dx1;
    if (delta === 0) return false; // parallel segments

    let s = (dx1 * (y21 - y11) + dy1 * (x11 - x21)) / delta;
    let t = (dx2 * (y11 - y21) + dy2 * (x21 - x11)) / -delta;

    return 0 <= s && s <= 1 && 0 <= t && t <= 1;
  }

  function point_segment_distance(px, py, x1, y1, x2, y2) {
    let dx = x2 - x1;
    let dy = y2 - y1;

    if (dx === 0 && dy === 0) {
      // the segment's just a point
      return Math.hypot(px - x1, py - y1);
    }

    // Calculate the t that minimizes the distance.
    let t = ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy);

    // See if this represents one of the segment's
    // end points or a point in the middle.
    if (t < 0) {
      dx = px - x1;
      dy = py - y1;
    } else if (t > 1) {
      dx = px - x2;
      dy = py - y2;
    } else {
      let near_x = x1 + t * dx;
      let near_y = y1 + t * dy;
      dx = px - near_x;
      dy = py - near_y;
    }
    return Math.hypot(dx, dy);
  }

  function areCollinear(vector1, vector2) {
    // Ensure vectors are of the same dimension
    if (vector1.length !== vector2.length) {
      return false;
    }

    // Find the ratio of the first non-zero pair of elements
    let ratio;
    for (let i = 0; i < vector1.length; i++) {
      if (vector1[i] !== 0 && vector2[i] !== 0) {
        ratio = vector1[i] / vector2[i];
        break;
      }
    }

    // Check if all corresponding elements are in the same ratio
    for (let i = 0; i < vector1.length; i++) {
      // Handle division by zero cases
      if (vector1[i] === 0 && vector2[i] !== 0 || vector1[i] !== 0 && vector2[i] === 0) {
        return false;
      }

      // Check the ratio
      if (vector1[i] !== 0 && vector2[i] !== 0) {
        if (vector1[i] / vector2[i] !== ratio) {
          return false;
        }
      }
    }

    return true;
  }

  /*  -----------------------  */
  /*  TODO modify lines below  */
  /*  -----------------------  */
  function distanceConstraint(agent_i, agent_j, desiredDistance) {
    const agentCentroidDist = distance(
      agent_i.px,
      agent_i.pz,
      agent_j.px,
      agent_j.pz
    );
    const agentDist = agentCentroidDist - desiredDistance;
    const dir_x = (agent_j.px - agent_i.px) / agentCentroidDist;
    const dir_z = (agent_j.pz - agent_i.pz) / agentCentroidDist;
    const agent_i_scaler =
      ((0.1 * agent_i.invmass) / (agent_i.invmass + agent_j.invmass)) *
      agentDist;
    const agent_j_scaler =
      ((0.1 * agent_j.invmass) / (agent_i.invmass + agent_j.invmass)) *
      agentDist;
    if (Math.abs(agentDist) > epsilon) {
      agent_i.px += agent_i_scaler * dir_x;
      agent_i.pz += agent_i_scaler * dir_z;
      agent_j.px += -agent_j_scaler * dir_x;
      agent_j.pz += -agent_j_scaler * dir_z;
    }
  }

  function collisionConstraint(agent_i, agent_j) {
    const agentCentroidDist = distance(
      agent_i.px,
      agent_i.pz,
      agent_j.px,
      agent_j.pz
    );
    const agentDist = agentCentroidDist - AGENTSIZE;
    const dir_x = (agent_j.px - agent_i.px) / agentCentroidDist;
    const dir_z = (agent_j.pz - agent_i.pz) / agentCentroidDist;
    const agent_i_scaler =
      (agent_i.invmass / (agent_i.invmass + agent_j.invmass)) * agentDist;
    const agent_j_scaler =
      (agent_j.invmass / (agent_i.invmass + agent_j.invmass)) * agentDist;
    if (agentDist < 0) {
      agent_i.px += agent_i_scaler * dir_x;
      agent_i.pz += agent_i_scaler * dir_z;
      agent_j.px += -agent_j_scaler * dir_x;
      agent_j.pz += -agent_j_scaler * dir_z;
    }
  }

  function agentVelocityPlanner() {
    sceneEntities.forEach(function (agent_i) {
      const distToGoal = distance(
        agent_i.x,
        agent_i.z,
        agent_i.goal_x,
        agent_i.goal_z
      );
      if (distToGoal > RADIUS) {
        const dir_x = (agent_i.goal_x - agent_i.x) / distToGoal;
        const dir_z = (agent_i.goal_z - agent_i.z) / distToGoal;
        agent_i.vx = agent_i.v_pref * dir_x;
        agent_i.vz = agent_i.v_pref * dir_z;
      }
      agent_i.vx = 0.9999 * agent_i.vx;
      agent_i.vz = 0.9999 * agent_i.vz;
    });
  }



  // function longRangeConstraint(agent_i, agent_j, agent_i_rad= AGENTSIZE, agent_j_rad = AGENTSIZE) {
  //   const agentCentroidDist = distance(agent_i.px, agent_i.pz,
  //       agent_j.px, agent_j.pz);
  //   const radius_init = agent_i_rad + agent_j_rad;
  //   const radius_sq_init = radius_init * radius_init;
  //   var radius_sq = radius_sq_init;
  //   const dv_i = 1.;  // 1./delta_t;
  //   let delta_correction_i = {"x":0, "y":0};
  //   let delta_correction_j= {"x":0, "y":0};
  //   if (agentCentroidDist < radius_init) {
  //     radius_sq = (radius_init - agentCentroidDist) * (radius_init - agentCentroidDist);
  //   }
  //   const v_x = (agent_i.px - agent_i.x) / timestep - (agent_j.px - agent_j.x) / timestep;
  //   const v_y = (agent_i.pz - agent_i.z) / timestep - (agent_j.pz - agent_j.z) / timestep;
  //   const x0 = agent_i.x - agent_j.x;
  //   const y0 = agent_i.z - agent_j.z;
  //   const v_sq = v_x * v_x + v_y * v_y;
  //   const x0_sq = x0 * x0;
  //   const y0_sq = y0 * y0;
  //   const x_sq = x0_sq + y0_sq;
  //   const a = v_sq;
  //   const b = -v_x * x0 - v_y * y0;   // b = -1 * v_.dot(x0_).  Have to check this.
  //   const b_sq = b * b;
  //   const c = x_sq - radius_sq;
  //   const d_sq = b_sq - a * c;
  //   const d = Math.sqrt(d_sq);
  //   const tao = (b - d) / a;
  //   let taoCalc = tao;
  //   // if (agent_i.index % 2 !== 0 && agent_j.index % 2 !== 0){
  //   //
  //   //   C_TAO0 = 6;
  //   // }else if (agent_i.index % 2 === 0 && agent_j.index % 2 === 0){
  //   //
  //   //   C_TAO0 = 80;
  //   // }else {
  //   //   return;
  //   // }
  //
  //   let lengthV;
  //   if (d_sq > 0.0 && Math.abs(a) > epsilon && tao > 0 && tao < C_TAU_MAX) {
  //     taoCalc=timestep* (1+  Math.floor(tao)/timestep);
  //     const c_tao = Math.exp(-tao * tao / C_TAO0);  //Math.abs(tao - C_TAO0);
  //     const tao_sq = c_tao * c_tao;
  //     const grad_x_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_x * taoCalc) - (x0 + (v_y * x0 * y0 + v_x * (radius_sq - y0_sq)) / d)));
  //     const grad_y_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_y * taoCalc) - (y0 + (v_x * x0 * y0 + v_y * (radius_sq - x0_sq)) / d)));
  //     const grad_x_j = -grad_x_i;
  //     const grad_y_j = -grad_y_i;
  //
  //
  //
  //
  //     const stiff = C_LONG_RANGE_STIFF * Math.exp(-tao * tao / C_TAO0);    //changed
  //     const s = stiff * tao_sq / (agent_i.invmass * (grad_y_i * grad_y_i + grad_x_i * grad_x_i) + agent_j.invmass * (grad_y_j * grad_y_j + grad_x_j * grad_x_j));     //changed
  //
  //
  //
  //     delta_correction_i = clamp2D(s * agent_i.invmass * grad_x_i,
  //         s * agent_i.invmass * grad_y_i,
  //         MAX_DELTA);
  //
  //     delta_correction_j = clamp2D(s * agent_j.invmass * grad_x_j,
  //         s * agent_j.invmass * grad_y_j,
  //         MAX_DELTA);
  //     agent_i.px += delta_correction_i.x;
  //     agent_i.pz += delta_correction_i.y;
  //     agent_j.px += delta_correction_j.x;
  //     agent_j.pz += delta_correction_j.y;
  //
  //     agent_i.grad[0] += delta_correction_i.x;
  //     agent_i.grad[1] += delta_correction_i.y;
  //     agent_j.grad[0] += delta_correction_j.x;
  //     agent_j.grad[1] += delta_correction_j.y;
  //
  //   }
  // }

  function longRangeConstraint(agent_i, agent_j, agent_i_rad = RADIUS, agent_j_rad = RADIUS) {
    const agentCentroidDist = distance(agent_i.px, agent_i.pz,
        agent_j.px, agent_j.pz);
    const radius_init = agent_i_rad + agent_j_rad // was AGENTSIZE
    const radius_sq_init = radius_init * radius_init;
    var radius_sq = radius_sq_init;
    const dv_i = 1.;  // 1./delta_t;
    let delta_correction_i = {"x":0, "y":0};
    let delta_correction_j= {"x":0, "y":0};
    if (agentCentroidDist < radius_init) {
      radius_sq = (radius_init - agentCentroidDist) * (radius_init - agentCentroidDist);
    }
    const v_x = (agent_i.px - agent_i.x) / timestep - (agent_j.px - agent_j.x) / timestep;
    const v_y = (agent_i.pz - agent_i.z) / timestep - (agent_j.pz - agent_j.z) / timestep;
    let x0 = agent_i.x - agent_j.x;
    let y0 = agent_i.z - agent_j.z;
    const v_sq = v_x * v_x + v_y * v_y;
    const y0_sq = y0 * y0;
    const x0_sq = x0 * x0;
    const x_sq = x0_sq + y0_sq;
    const a = v_sq;
    const b = -v_x * x0 - v_y * y0;   // b = -1 * v_.dot(x0_).  Have to check this.
    const b_sq = b * b;
    const c = x_sq - radius_sq;
    const d_sq = b_sq - a * c;
    const d = Math.sqrt(d_sq);
    const tao = (b - d) / a;
    // if (agent_i.index % 2 !== 0 && agent_j.index % 2 !== 0){
    //
    //   C_TAO0 = 6;
    // }else if (agent_i.index % 2 === 0 && agent_j.index % 2 === 0){
    //
    //   C_TAO0 = 80;
    // }else {
    //   return;
    // }
    let lengthV;
    if (d_sq > 0.0 && tao > 0 && tao < C_TAU_MAX) {
      let v_minus_x_diff_theta = b / ( Math.sqrt(v_sq) *  Math.sqrt(x0_sq))
      let v_diff_theta = Math.atan2(v_y, v_x); //v_y/v_x
      let x_diff_theta = Math.atan2(y0, x0);//y0/x0

      if(Math.abs(v_diff_theta-x_diff_theta) > 0.999*Math.PI &&
          Math.abs(v_diff_theta-x_diff_theta) < 1.001*Math.PI
      )
      {
        let r_xdiff = Math.sqrt(x_sq);
        let x0_mod = r_xdiff * Math.cos(x_diff_theta)
        let y0_mod = r_xdiff * Math.sin(x_diff_theta)
        x0 = r_xdiff * Math.cos(x_diff_theta + Math.PI*0.001);
        y0 = r_xdiff * Math.sin(x_diff_theta + Math.PI*0.001);
        console.log("v-x angel " + (v_diff_theta-x_diff_theta) + " xangle " + x_diff_theta )
        console.log("changed: x0,y0  (" + x0 + "," + y0 + ") new:" + x0_mod + " " +y0_mod);
      }
      const c_tao = Math.exp(-tao * tao / C_TAO0);  //Math.abs(tao - C_TAO0);
      const tao_sq = c_tao * c_tao;
      const grad_x_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_x * tao) - (x0 + (v_y * x0 * y0 + v_x * (radius_sq - y0_sq)) / d)));
      const grad_y_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_y * tao) - (y0 + (v_x * x0 * y0 + v_y * (radius_sq - x0_sq)) / d)));
      const grad_x_j = -grad_x_i;
      const grad_y_j = -grad_y_i;
      const stiff = C_LONG_RANGE_STIFF * Math.exp(-tao * tao / C_TAO0);    //changed
      //console.log("tao " + tao + " stiff " + stiff );
      const s = stiff * tao_sq / (agent_i.invmass * (grad_y_i * grad_y_i + grad_x_i * grad_x_i) + agent_j.invmass * (grad_y_j * grad_y_j + grad_x_j * grad_x_j));     //changed

      lengthV = Math.sqrt(s * agent_i.invmass * grad_x_i * s * agent_i.invmass * grad_x_i
          + s * agent_i.invmass * grad_y_i * s * agent_i.invmass * grad_y_i);

      delta_correction_i = clamp2D(s * agent_i.invmass * grad_x_i,
          s * agent_i.invmass * grad_y_i,
          MAX_DELTA);

      delta_correction_j =  clamp2D(s * agent_j.invmass * grad_x_j,
          s * agent_j.invmass * grad_y_j,
          MAX_DELTA);
      agent_i.px +=  delta_correction_i.x;
      agent_i.pz += 44* delta_correction_i.y;
      agent_j.px += delta_correction_j.x;
      agent_j.pz += delta_correction_j.y;

      agent_i.grad[0] += delta_correction_i.x;
      agent_i.grad[1] += delta_correction_i.y;
      agent_j.grad[0] += delta_correction_j.x;
      agent_j.grad[1] += delta_correction_j.y;

      // for utilities
      agent_i.grad.x += grad_x_i;
      agent_i.grad.z += grad_y_i;
      agent_j.grad.x += grad_x_j;
      agent_j.grad.z += grad_y_j;

    }
  }


  let C_TAU_MAX = 20;
  // const C_MAX_ACCELERATION = 0.01;
  let C_TAO0 = 20; //
  const C_LONG_RANGE_STIFF = 1.0;

  function clamp2D(vx,vy, maxValue) {
    const lengthV = Math.sqrt(vx * vx + vy * vy);
    if (lengthV > maxValue) {
      const mult = (maxValue / lengthV);
      vx *= mult;
      vy *= mult;
    }
    return {"x":vx, "y":vy}
  }

  function longRangeConstraintCapsuleV0(best_i, best_j,
                                      p_best_i, p_best_j,
                                      theta_i, theta_j,
                                      agent_i, agent_j,
                                      entity_i, entity_j,
                                      i = -1, j = -1) {


    const agentCentroidDist = distance(p_best_i.x, p_best_i.z, p_best_j.x, p_best_j.z);

    const radius_init = 2 * AGENTSIZE;
    const radius_sq_init = radius_init * radius_init;
    let radius_sq = radius_sq_init;
    const dv_i = 1.;  // 1./delta_t;
    let delta_correction_i = {"x":0, "y":0};
    let delta_correction_j= {"x":0, "y":0};
    if (agentCentroidDist < radius_init) {
      radius_sq = (radius_init - agentCentroidDist) * (radius_init - agentCentroidDist);
    }

    let v_x = (p_best_i.x - best_i.x) / timestep - (p_best_j.x - best_j.x) / timestep;
    let v_y = (p_best_i.z - best_i.z) / timestep - (p_best_j.z - best_j.z) / timestep;
    let x0 = best_i.x - best_j.x;
    let y0 = best_i.z - best_j.z;



    const v_sq = v_x * v_x + v_y * v_y;
    const x0_sq = x0 * x0;
    const y0_sq = y0 * y0;
    const x_sq = x0_sq + y0_sq;
    const a = v_sq;
    const b = -v_x * x0 - v_y * y0;   // b = -1 * v_.dot(x0_).  Have to check this.
    const b_sq = b * b;
    const c = x_sq - radius_sq;
    const d_sq = b_sq - a * c;
    const d = Math.sqrt(d_sq);
    const tao = (b - d) / a;
    // console.log("ttc in long range paper: " + tao);


    let dx = 0;
    let dz = 0;

    if(C_TAU_MAX <= tao < 5 + C_TAU_MAX ){

    }



    let grad_x_i;
    let grad_y_i;
    let grad_x_j;
    let grad_y_j;
    let s;

    if (d_sq > 0.0 && Math.abs(a) > epsilon && tao > 0 && tao < C_TAU_MAX) {
      const c_tao = Math.exp(-tao * tao / C_TAO0);  //Math.abs(tao - C_TAO0);
      const tao_sq = c_tao * c_tao;

      grad_x_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_x * tao) - (x0 + (v_y * x0 * y0 + v_x * (radius_sq - y0_sq)) / d)));
      grad_y_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_y * tao) - (y0 + (v_x * x0 * y0 + v_y * (radius_sq - x0_sq)) / d)));
      grad_x_j = -grad_x_i;
      grad_y_j = -grad_y_i;

      // special case
      console.log(i + "<==>" + j)




      // if facing direction on the same line AND the best points are exactly facing with each other
      // adding gradient value



      const stiff = C_LONG_RANGE_STIFF * Math.exp(-tao * tao / C_TAO0);    //changed
      s = stiff * tao_sq / (0.5 * (grad_y_i * grad_y_i + grad_x_i * grad_x_i) + 0.5 * (grad_y_j * grad_y_j + grad_x_j * grad_x_j));     //changed
      // console.log()


      delta_correction_i = clamp2D(s * 0.5 * grad_x_i,
          s * 0.5 * grad_y_i,
          MAX_DELTA);

      delta_correction_j = clamp2D(s * 0.5 * grad_x_j,
          s * 0.5 * grad_y_j,
          MAX_DELTA);

      // console.log("Long Range active");

    }else {
      grad_x_i = 0;
      grad_y_i = 0;
      grad_x_j = 0;
      grad_y_j = 0;
      s=0;
    }


    // return tao;
    return [
      delta_correction_i,
      delta_correction_j,
      [grad_x_i, grad_y_i],
      [grad_x_j, grad_y_j],
      s
    ];
  }


  function longRangeConstraintCapsule(best_i, best_j,
                                      p_best_i, p_best_j,
                                      theta_i, theta_j,
                                      agent_i, agent_j,
                                      entity_i, entity_j,
                                      i = -1, j = -1,
                                      correctX = 0, correctZ = 0
                                      ) {


    const agentCentroidDist = distance(p_best_i.x, p_best_i.z, p_best_j.x, p_best_j.z);

    const radius_init = 2 * AGENTSIZE;
    const radius_sq_init = radius_init * radius_init;
    let radius_sq = radius_sq_init;
    const dv_i = 1.;  // 1./delta_t;
    let delta_correction_i = {"x":0, "y":0};
    let delta_correction_j= {"x":0, "y":0};
    if (agentCentroidDist < radius_init) {
      radius_sq = (radius_init - agentCentroidDist) * (radius_init - agentCentroidDist);
    }

    let v_x = (p_best_i.x - best_i.x) / timestep - (p_best_j.x - best_j.x) / timestep;
    let v_y = (p_best_i.z - best_i.z) / timestep - (p_best_j.z - best_j.z) / timestep;
    let x0 = best_i.x - best_j.x;
    let y0 = best_i.z - best_j.z;


    // rotation (facing angle) difference
    // best points distance difference

    const v_sq = v_x * v_x + v_y * v_y;
    const x0_sq = x0 * x0;
    const y0_sq = y0 * y0;
    const x_sq = x0_sq + y0_sq;
    const a = v_sq;
    const b = -v_x * x0 - v_y * y0;   // b = -1 * v_.dot(x0_).  Have to check this.
    const b_sq = b * b;
    const c = x_sq - radius_sq;
    const d_sq = b_sq - a * c;
    const d = Math.sqrt(d_sq);
    const tao = (b - d) / a;
    // console.log("ttc in long range paper: " + tao);
    let taoCalc = tao;
    let grad_x_i;
    let grad_y_i;
    let grad_x_j;
    let grad_y_j;
    let s;


    if (d_sq > 0.0 && Math.abs(a) > epsilon && tao > 0 && tao < C_TAU_MAX) {

      customParams.status[i][j] = true
      customParams.status[j][i] = true

      let v_diff_theta = Math.atan2(v_y, v_x); //v_y/v_x
      let x_diff_theta = Math.atan2(y0, x0);//y0/x0
      
      if(Math.abs(v_diff_theta-x_diff_theta) > 0.99*Math.PI && 
      Math.abs(v_diff_theta-x_diff_theta) < 1.01*Math.PI 
      )
      {
          let r_xdiff = Math.sqrt(x_sq);
          let x0_mod = r_xdiff * Math.cos(x_diff_theta)
          let y0_mod = r_xdiff * Math.sin(x_diff_theta)
          x0 = r_xdiff * Math.cos(x_diff_theta + Math.PI*0.01);
          y0 = r_xdiff * Math.sin(x_diff_theta + Math.PI*0.01);
          console.log("v-x angel " + (v_diff_theta-x_diff_theta) + " xangle " + x_diff_theta )
          console.log("changed: x0,y0  (" + x0 + "," + y0 + ") new:" + x0_mod + " " +y0_mod);
      }

      taoCalc = (Math.floor(tao)/timestep) * timestep;
      const c_tao = Math.exp(-taoCalc * taoCalc / C_TAO0);  //Math.abs(tao - C_TAO0);
      const tao_sq = c_tao * c_tao;

      taoCalc+=timestep;
      grad_x_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_x * taoCalc) - (x0 + (v_y * x0 * y0 + v_x * (radius_sq - y0_sq)) / d))) ;
      grad_y_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_y * taoCalc) - (y0 + (v_x * x0 * y0 + v_y * (radius_sq - x0_sq)) / d)));
      grad_x_j = -grad_x_i;
      grad_y_j = -grad_y_i;

      // special case
      console.log(i + "<==>" + j)

      // if facing direction on the same line AND the best points are exactly facing with each other
      // adding gradient value

      const stiff = C_LONG_RANGE_STIFF * Math.exp(-tao * tao / C_TAO0);    //changed
      s = stiff * tao_sq / (0.5 * (grad_y_i * grad_y_i + grad_x_i * grad_x_i) + 0.5 * (grad_y_j * grad_y_j + grad_x_j * grad_x_j));     //changed

      delta_correction_i = clamp2D(s * 0.5 * grad_x_i + s * 0.5 * -correctX,
          s * 0.5 * grad_y_i + s * 0.5 * -correctZ,
          MAX_DELTA);

      delta_correction_j = clamp2D(s * 0.5 * grad_x_j + s * 0.5 * correctX,
          s * 0.5 * grad_y_j + s * 0.5 * correctZ,
          MAX_DELTA) ;

      // console.log("Long Range active");

    }else {
      grad_x_i = 0;
      grad_y_i = 0;
      grad_x_j = 0;
      grad_y_j = 0;
      s=0;
    }


    // return tao;
    return [
        delta_correction_i,
        delta_correction_j,
        [grad_x_i, grad_y_i],
        [grad_x_j, grad_y_j],
        s
    ];
  }

  function getBestPointWithWall(xi, zi, wall){

    const iCoords = rotateLineSegment(
        xi,
        zi + agentLength + RADIUS,
        xi,
        zi - agentLength - RADIUS,
        sceneEntities[i].agent.rotation.z
    );



    // Agent A
    const a = {
      tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
      base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
      radius: RADIUS,
    };
    // Wall B
    const b = {
      tip: new THREE.Vector3(wall.tip.x, 0, wall.tip.z),
      base: new THREE.Vector3(wall.base.x, 0, wall.base.z),
      radius: RADIUS,
    };


    // capsule A:
    const a_Normal = a.tip.clone().sub(a.base.clone()).normalize();
    const a_LineEndOffset = a_Normal.clone().multiplyScalar(a.radius);
    const a_A = a.base.clone().add(a_LineEndOffset);
    const a_B = a.tip.clone().sub(a_LineEndOffset);

    // capsule B:
    const b_Normal = b.tip.clone().sub(b.base.clone()).normalize();
    const b_LineEndOffset = b_Normal.clone().multiplyScalar(b.radius);
    const b_A = b.base.clone().add(b_LineEndOffset);
    const b_B = b.tip.clone().sub(b_LineEndOffset);

    // vectors between line endpoints:
    const v0 = b_A.clone().sub(a_A);
    const v1 = b_B.clone().sub(a_A);
    const v2 = b_A.clone().sub(a_B);
    const v3 = b_B.clone().sub(a_B);

    // squared distances:
    const d0 = v0.clone().dot(v0);
    const d1 = v1.clone().dot(v1);
    const d2 = v2.clone().dot(v2);
    const d3 = v3.clone().dot(v3);

    // select best potential endpoint on capsule A:
    let bestA;
    if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1) {
      bestA = a_B;
    } else {
      bestA = a_A;
    }

    // select point on capsule B line segment nearest to best potential endpoint on A capsule:
    const bestB = ClosestPointOnLineSegment(b_A, b_B, bestA);

    // now do the same for capsule A segment:
    bestA = ClosestPointOnLineSegment(a_A, a_B, bestB);

    return [bestA, bestB, a, b]
  }



  function getBestPoint(xi, zi, xj, zj){

    const iCoords = rotateLineSegment(
        xi,
        zi + agentLength + RADIUS,
        xi,
        zi - agentLength - RADIUS,
        sceneEntities[i].agent.rotation.z
    );

    const jCoords = rotateLineSegment(
        xj,
        zj + agentLength + RADIUS,
        xj,
        zj - agentLength - RADIUS,
        sceneEntities[j].agent.rotation.z
    );


    // Agent A
    const a = {
      tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
      base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
      radius: RADIUS,
      real_tip: null,
      real_base: null
    };
    // Agent B
    const b = {
      tip: new THREE.Vector3(jCoords[0], 0, jCoords[1]),
      base: new THREE.Vector3(jCoords[2], 0, jCoords[3]),
      radius: RADIUS,
      real_tip: null,
      real_base: null
    };


    // capsule A:
    const a_Normal = a.tip.clone().sub(a.base.clone()).normalize();
    const a_LineEndOffset = a_Normal.clone().multiplyScalar(a.radius);
    const a_A = a.base.clone().add(a_LineEndOffset);
    const a_B = a.tip.clone().sub(a_LineEndOffset);
    a.real_tip = a_B;
    a.real_base = a_A;

    // capsule B:
    const b_Normal = b.tip.clone().sub(b.base.clone()).normalize();
    const b_LineEndOffset = b_Normal.clone().multiplyScalar(b.radius);
    const b_A = b.base.clone().add(b_LineEndOffset);
    const b_B = b.tip.clone().sub(b_LineEndOffset);
    b.real_tip = b_B;
    b.real_base = b_A;

    // vectors between line endpoints:
    const v0 = b_A.clone().sub(a_A);
    const v1 = b_B.clone().sub(a_A);
    const v2 = b_A.clone().sub(a_B);
    const v3 = b_B.clone().sub(a_B);

    // squared distances:
    const d0 = v0.clone().dot(v0);
    const d1 = v1.clone().dot(v1);
    const d2 = v2.clone().dot(v2);
    const d3 = v3.clone().dot(v3);

    // select best potential endpoint on capsule A:
    let bestA;
    if (d2 < d0 || d2 < d1 || d3 < d0 || d3 < d1) {
      bestA = a_B;
    } else {
      bestA = a_A;
    }

    // select point on capsule B line segment nearest to best potential endpoint on A capsule:
    const bestB = ClosestPointOnLineSegment(b_A, b_B, bestA);

    // now do the same for capsule A segment:
    bestA = ClosestPointOnLineSegment(a_A, a_B, bestB);

    return [bestA, bestB, a, b]
  }



  function checkAgentAlignment(agentA, agentB){

    let projectBaseB = ClosestPointOnLineSegment(agentA.base, agentA.tip, agentB.base);
    let projectTipB = ClosestPointOnLineSegment(agentA.base, agentA.tip, agentB.tip);

    // calculate distance difference case 1
    let diff1 = projectBaseB.clone().sub(agentA.base)
    let diff2 = projectTipB.clone().sub(agentA.tip)

    // calculate distance difference case 2
    let diff3 = projectTipB.clone().sub(agentA.base)
    let diff4 = projectBaseB.clone().sub(agentA.tip)

    // return (diff1.length() < 0.5 && diff2.length() < 0.5) || (diff3.length() < 0.5 && diff4.length() < 0.5);
    return (diff1.length() < 0.1 && diff2.length() < 0.1) || (diff3.length() < 0.1 && diff4.length() < 0.1);

  }


  function findLargestAngleVectors(capsule1, capsule2) {
    // capsule1 and capsule2 are objects with {tip: THREE.Vector3, base: THREE.Vector3}
    // note that the capsule2 is the expanded capsule

    // Calculate vectors between tips and bases
    let tipToTip = capsule2.tip.clone().sub(capsule1.tip);
    let baseToBase = capsule2.base.clone().sub(capsule1.base);
    let tipToBase = capsule2.base.clone().sub(capsule1.tip);
    let baseToTip = capsule2.tip.clone().sub(capsule1.base);

    // Calculate angles between dir1 and other vectors
    let angle1 = baseToBase.angleTo(tipToTip);
    let angle2 = baseToBase.angleTo(tipToBase);
    let angle3 = baseToBase.angleTo(baseToTip);

    let angle4 = tipToTip.angleTo(tipToBase);
    let angle5 = tipToTip.angleTo(baseToTip);

    let angle6 = tipToBase.angleTo(baseToTip);


    // Identify the pair of vectors with the largest angle
    let maxAngle = Math.max(angle1, angle2, angle3, angle4, angle5, angle6);
    let vectorsWithLargestAngle;

    if (maxAngle === angle1) {
      vectorsWithLargestAngle = {vector1: baseToBase, vector2: tipToTip};
    } else if (maxAngle === angle2) {
      vectorsWithLargestAngle = {vector1: baseToBase, vector2: tipToBase};
    } else if (maxAngle === angle3) {
      vectorsWithLargestAngle = {vector1: baseToBase, vector2: baseToTip};
    } else if (maxAngle === angle4) {
      vectorsWithLargestAngle = {vector1: tipToTip, vector2: tipToBase};
    }else if (maxAngle === angle5) {
      vectorsWithLargestAngle = {vector1: tipToTip, vector2: baseToTip};
    }else {
      vectorsWithLargestAngle = {vector1: tipToBase, vector2: baseToTip};
    }

    return vectorsWithLargestAngle;
  }

  function formNewCapsules(xi, zi, xj, zj, rad_i = RADIUS, rad_j = RADIUS){
    const iCoords = rotateLineSegment(
        xi,
        zi + agentLength + RADIUS,
        xi,
        zi - agentLength - RADIUS,
        sceneEntities[i].agent.rotation.z
    );

    const jCoords = rotateLineSegment(
        xj,
        zj + agentLength + RADIUS,
        xj,
        zj - agentLength - RADIUS,
        sceneEntities[j].agent.rotation.z
    );


    // Agent A
    const a = {
      center: new THREE.Vector3(xi, 0, zi),
      tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
      base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
      radius: RADIUS,
      real_tip: null,
      real_base: null
    };
    // Agent B
    const b = {
      center: new THREE.Vector3(xj, 0, zj),
      tip: new THREE.Vector3(jCoords[0], 0, jCoords[1]),
      base: new THREE.Vector3(jCoords[2], 0, jCoords[3]),
      radius: RADIUS,
      real_tip: null,
      real_base: null
    };

    return [a, b]
  }

  function projectPointOntoVector(point, vector) {
    // Clone the vector to avoid mutating the original
    let direction = vector.clone().normalize(); // Ensure the vector is a unit vector
    let pointVector = point.clone(); // Clone to avoid mutating the original point

    // Calculate the dot product between pointVector and direction
    let dotProduct = pointVector.dot(direction);

    // Calculate the projection of pointVector onto the direction vector
    // proj = (dotProduct / direction.lengthSq()) * direction
    return direction.multiplyScalar(dotProduct / direction.lengthSq()); // This is the projected point as a Vector3
  }


  function getBestBDPrime(base, tip, velocity){
    let B_candidate1 =  base.clone().add(velocity.clone().multiplyScalar(RADIUS));
    let B_candidate2 =  base.clone().sub(velocity.clone().multiplyScalar(RADIUS));

    let T_candidate1 =  tip.clone().add(velocity.clone().multiplyScalar(RADIUS));
    let T_candidate2 =  tip.clone().sub(velocity.clone().multiplyScalar(RADIUS));

    let dist1 = B_candidate1.distanceTo(T_candidate1);
    let dist2 = B_candidate1.distanceTo(T_candidate2);
    let dist3 = B_candidate2.distanceTo(T_candidate1);
    let dist4 = B_candidate2.distanceTo(T_candidate2);

    let maxDist = Math.max(dist1, dist2, dist3, dist4);

    if (maxDist === dist1){
      return [B_candidate1, T_candidate1, maxDist]

    }else if(maxDist === dist2){
      return [B_candidate1, T_candidate2, maxDist]

    }else if(maxDist === dist3){
      return [B_candidate2, T_candidate1, maxDist]

    }else {
      return [B_candidate2, T_candidate2, maxDist]

    }
  }

  function getPerpendicular2DVector(vec) {
    // Given a vector {x, z}, returns a vector that is perpendicular to it
    return new THREE.Vector3(-vec.z, 0, vec.x);
  }

  function checkVectorDirection(vectorA, vectorB) {
    // Assuming vectorA and vectorB are instances of THREE.Vector3
    const dotProduct = vectorA.dot(vectorB);

    if (dotProduct > 0) {
      return 1;
    } else if (dotProduct < 0) {
      return -1;
    } else {
      return 0;
    }
  }

  function determineHeadTail(velPrime, basePrime, tipPrime, centerPrime, another_centerPrime){

    let head, tail;

    let flag_i = checkVectorDirection(velPrime, basePrime.clone().sub(centerPrime))
    if(flag_i === 1){
      head = basePrime;
      tail = tipPrime;
    }else if(flag_i === -1){
      head = tipPrime;
      tail = basePrime;
    }else {

      let dist1 = basePrime.distanceTo(another_centerPrime);
      let dist2 = tipPrime.distanceTo(another_centerPrime);

      let minDist = Math.min(dist1, dist2);
      if(minDist === dist1){
        head = basePrime;
        tail = tipPrime;

      }else {
        head = tipPrime;
        tail = basePrime;

      }
    }

    return [head, tail];
  }

  function timeToCollideVector3(P1, V1, P2, V2) {
    // Assuming linear motion and choosing one component to represent the line direction, e.g., the x-component
    let deltaP = P2.x - P1.x; // Difference in positions
    let deltaV = V1.x - V2.x; // Difference in velocities

    // Check if velocities are equal in the chosen direction
    if (deltaV === 0) {
      // If starting at the same position, they collide immediately, otherwise never
      return deltaP === 0 ? 0 : null;
    } else {
      // Calculate time to collision based on the chosen component
      // If time is negative, they are moving apart or have already passed each other
      return deltaP / deltaV;
    }
  }



  /*  -----------------------  */

  agentVelocityPlanner();

  sceneEntities.forEach(function (item) {
    item.px = item.x + timestep * item.vx;
    item.pz = item.z + timestep * item.vz;
    item.py = item.y + timestep * item.vy;
  });






  let pbdIters = 0;
  let isColliding;
  let agent_a,
    agent_b,
    desDistance,
    i,
    j,
    idx = 0;

  while (pbdIters < ITERNUM) {


    // clean previous accumulated gradient
    i = 0;
    while (i < sceneEntities.length) {
      j = i + 1;
      while (j < sceneEntities.length) {

        sceneEntities[i].grad.x = 0;
        sceneEntities[i].grad.z = 0;
        sceneEntities[j].grad.x = 0;
        sceneEntities[j].grad.z = 0;

        sceneEntities[i].grad.dx = 0;
        sceneEntities[i].grad.dz = 0;
        sceneEntities[j].grad.dx = 0;
        sceneEntities[j].grad.dz = 0;

        customParams.status[i][j] = false;
        customParams.status[j][i] = false;

        j += 1;
      }
      i += 1;
    }

    // wall collision (based on short range)
    i=0;
    while(i<sceneEntities.length)
    {
      j=0;
      while(j<customParams.wallData.length)
      {
        let [p_bestA, w_bestB, p_agent_i,p_agent_j] = getBestPointWithWall(sceneEntities[i].px, sceneEntities[i].pz, customParams.wallData[j]);

        let penetration_normal = p_bestA.clone().sub(w_bestB);
        const len = penetration_normal.length();
        penetration_normal.divideScalar(len); // normalize
        const penetration_depth = sceneEntities[i].radius + Math.sqrt(2) * sceneEntities[i].radius - len;
        const intersects = penetration_depth > 0;
        if (intersects) {
          sceneEntities[i].colliding = true;

          sceneEntities[i].px += penetration_normal.x * 1 * penetration_depth;
          sceneEntities[i].pz += penetration_normal.z * 1 * penetration_depth;

        }

        j+=1;
      }
      i+=1
    }


    i = 0;
    while (i < sceneEntities.length) {
      j = i + 1;
      while (j < sceneEntities.length) {

        // get capsule shape
        let [capsule_i, capsule_j] = formNewCapsules(sceneEntities[i].x, sceneEntities[i].z, sceneEntities[j].x, sceneEntities[j].z);

        let vi = new THREE.Vector3(sceneEntities[i].vx, 0, sceneEntities[i].vz);
        let vj = new THREE.Vector3(sceneEntities[j].vx, 0, sceneEntities[j].vz);

        // we need to use these 4 axes to find if they are separate or not
        let axes = [vi, vj, getPerpendicular2DVector(vi), getPerpendicular2DVector(vj)];
        let ttc = [];
        for (let k = 0; k<axes.length;k++){
          // check current axis
          let axis = axes[k];
          // projected tip and base of j on to velocity of i
          let base_j_prime =  projectPointOntoVector(capsule_j.base, axis)
          let tip_j_prime =  projectPointOntoVector(capsule_j.tip, axis)
          let center_j_prime = projectPointOntoVector(capsule_j.center, axis)
          let [base_j_, tip_j_, dist_j] = getBestBDPrime(base_j_prime, tip_j_prime, axis);

          // // projected tip and base of i on to velocity of j
          let base_i_prime =  projectPointOntoVector(capsule_i.base, axis)
          let tip_i_prime =  projectPointOntoVector(capsule_i.tip, axis)
          let center_i_prime = projectPointOntoVector(capsule_i.center, axis)
          let [base_i_, tip_i_, dist_i] = getBestBDPrime(base_i_prime, tip_i_prime, axis);

          let pVel_i, pVel_j;

          pVel_i = axis.clone().multiplyScalar(vi.dot(axis) / axis.lengthSq());
          pVel_j = axis.clone().multiplyScalar(vj.dot(axis) / axis.lengthSq());


          let [head_i, tail_i] = determineHeadTail(pVel_i, base_i_, tip_i_, center_i_prime, center_j_prime);
          let [head_j, tail_j] = determineHeadTail(pVel_j, base_j_, tip_j_, center_j_prime, center_i_prime);

          let t = timeToCollideVector3(head_i, pVel_i, head_j, pVel_j)
          ttc.push(t);

        }

        console.log(ttc);



        // longRangeConstraint(sceneEntities[i], sceneEntities[j], agentLength + RADIUS, agentLength + RADIUS )
        // sceneEntities[i].sphere = bestA;
        // sceneEntities[j].sphere = bestB;
        //
        // sceneEntities[i].px += delta_correction_i.x;
        // sceneEntities[i].pz += delta_correction_i.y;
        // sceneEntities[j].px += delta_correction_j.x;
        // sceneEntities[j].pz += delta_correction_j.y;
        //
        //
        // // for utilities
        // sceneEntities[i].grad.x += grad_i[0];
        // sceneEntities[i].grad.z += grad_i[1];
        // sceneEntities[j].grad.x += grad_j[0];
        // sceneEntities[j].grad.z += grad_j[1];
        //
        // sceneEntities[i].grad.s = s;
        // sceneEntities[j].grad.s = s;
        //
        // sceneEntities[i].grad.dx += delta_correction_i.x;
        // sceneEntities[i].grad.dz += delta_correction_i.y;
        // sceneEntities[j].grad.dx += delta_correction_j.x;
        // sceneEntities[j].grad.dz += delta_correction_j.y;
        //
        // customParams.best[i][j] = [bestA, bestB]
        // customParams.best[j][i] = [bestB, bestA]


        j += 1;
      }
      i += 1;
    }





    pbdIters += 1;
  }


  sceneEntities.forEach(function (item) {

    const dx = item.px - item.x;
    const dz = item.pz - item.z;
    item.agent.rotation.z = Math.atan2(dz, dx);

    item.vx = (item.px - item.x) / timestep;
    item.vz = (item.pz - item.z) / timestep;
    item.vy = (item.py - item.y) / timestep;
    item.x = item.px;
    item.z = item.pz;
    item.y = item.py;
  });
}

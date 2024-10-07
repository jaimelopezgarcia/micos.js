import * as math_utils from "./math_utils.js";
import * as math from 'mathjs'; 
import { getConstraintsRigid } from "./solver.js";
import {translate, rotate} from "./math_utils.js";

function compassWalkerState(){
    //Simple gated compas walker, basically 3 masses with 2 distance constraints
    // when one of the masses contacts the ground, its dynamics are the same as an inverted double pendulum
    //if springActuator is not null, it must be a SpringActuator object, and we create a cos(t) contraction function as in the example

    //lets add two extra masses as rigid knees, just to attach the springActuator

    let initialAngleDeg = 43;
    let initialAngle = initialAngleDeg*Math.PI/180;
    let uVec = [Math.cos(initialAngle), -Math.sin(initialAngle)];
    let legLength = 2.7;
    let xs = [[0,0],
              [0,legLength],
              [legLength*uVec[0],legLength*uVec[1]+legLength],
            ];
    
    
    
    //lets rotate the xs 20º around the first particle
    
    xs = micos.math_utils.rotate(xs,4,xs[0]);
    // lets translate the xs to the left
    xs = micos.math_utils.translate(xs,[-0.5,-0.3]);
    // now lets add the distance constraints between [[0,1],[0,2],[1,2],[2,3],[2,4],[3,4]]
    
    let constraints_distance = [[0,1,legLength],
                                [1,2,legLength]
                            ];
    
    let masses = [3,3,3];
    masses = masses.map((m) => m*1);
    let ground = [[-4,-1],[-4,-0.5],[4,-0.5],[4,-1]];
    ground = micos.math_utils.rotateCOM(ground,4);
    let dampingCoef = 0.000;
    //lets translate xs and ground to the left
    xs = micos.math_utils.translate(xs,[-2.5,-0.3]);
    ground = micos.math_utils.translate(ground,[-2.5,-0.3]);
    let polygons = [ground];

    let STATE = {
        "xs": xs,
        "vs": xs.map((x) => {return [0,0]}),
        "masses": masses,
        "constraints_distance": constraints_distance,
        "constraints_pin": [],
        "polygons": polygons,
        "t":0,
        "gravity": [1,[0,-1]],
        "springs": [],//particle1, particle2, stiffness, restLength
        "damping": xs.map((x,i) => [[i],dampingCoef]),

    };

    return STATE



}



function compassWalkerState2(groundAngle= 4.2){
    //Simple gated compas walker, basically 3 masses with 2 distance constraints
    // when one of the masses contacts the ground, its dynamics are the same as an inverted double pendulum
    //if springActuator is not null, it must be a SpringActuator object, and we create a cos(t) contraction function as in the example

    //lets add two extra masses as rigid knees, just to attach the springActuator

    let initialAngleDeg = 40;
    let initialAngle = initialAngleDeg*Math.PI/180;
    let uVec = [Math.cos(initialAngle), -Math.sin(initialAngle)];
    let legLength = 0.5;
    let xs = [[0,0],
              [0,legLength],
              [0,2*legLength],//head
              [legLength*uVec[0],legLength*uVec[1]+2*legLength],
              [2*legLength*uVec[0],2*legLength*uVec[1]+2*legLength]
            ];

    

    //lets rotate the xs 20º around the first particle

    xs = micos.math_utils.rotate(xs,4,xs[0]);
    // lets translate the xs to the left
    xs = micos.math_utils.translate(xs,[-0.5,-0.3]);
    // now lets add the distance constraints between [[0,1],[0,2],[1,2],[2,3],[2,4],[3,4]]

    let constraints_distance = [[0,1,legLength],
                                [0,2,2*legLength],
                                [1,2,legLength],
                                [2,3,legLength],
                                [2,4,2*legLength],
                                [3,4,legLength]];

    let masses = [3,1,3,1,3];


    // lets create a wide ground polygon that spans from -4 to 4 in x and -1 to -0.5 in y, and rotate it 20º

    // lets create the polygon in clockwise order
    let ground = [[-4,-1],[-4,-0.5],[4,-0.5],[4,-1]];
    ground = micos.math_utils.rotateCOM(ground, groundAngle);
    let polygons = [ground];
    let dampingCoef = 0.0001;

    let STATE = {
        "xs": xs,
        "vs": xs.map((x) => {return [0,0]}),
        "masses": masses,
        "constraints_distance": constraints_distance,
        "constraints_pin": [],
        "polygons": polygons,
        "t":0,
        "gravity": [0.5,[0,-1]],
        "springs": [],//particle1, particle2, stiffness, restLength
        "damping": xs.map((x,i) => [[i],dampingCoef]),

    };

    return STATE



}


function rigidDoublePendulum(){
    //lets make a double pendulum  but with a square replacing the masses

    let xs = [[0,0],[1,0],[1,-1],[0,-1]];
    let pinPoint = [0,1];
    let constraints_distance = getConstraintsRigid(xs);
    //lets add a triangle attached to  the 3 mass of the square

    let triangle = [[1,-2],[1,-3],[0,-2.5]];
    let extra_constraints_distance = [[2,4,1],[4,5,1],[5,6,1],[6,4,1]];
    constraints_distance = constraints_distance.concat(extra_constraints_distance);
    xs = xs.concat(triangle);


    //lets rotate xs
    xs = micos.math_utils.rotate(xs,130,xs[0]);
    xs = micos.math_utils.translate(xs,[0,1]);
    let masses = [1,1,1,1,1,1,1];
    let constraints_pin = [[0,pinPoint,1]];
    let gravity = [1,[0,-1]];
    let time = 0;
    let dampingCoef = 0.01;
    let initialState = {
        xs: xs,
        vs: [[0,0],[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]],
        masses: masses,
        constraints_distance: constraints_distance,
        constraints_pin: constraints_pin,
        gravity: gravity,
        time: time,
        damping : xs.map((x,i) => [[i],dampingCoef]),
        }

    return initialState;



}




function bungeeJumperState(){
    //a couple of masses for feet knees, a big hip and torso a head and two arms
    // the arms are pinned to the torso, and a spring is attached to the feet ( upside down)

    let feetL = [-0.5,0];
    let kneeL = [-0.5,-1];
    let feetR = [0.5,0];
    let kneeR = [0.5,-1];
    let hip = [0,-2];
    let torso = [0,-3];
    let head = [0,-4];
    let armL = [-1,-3];
    let armR = [1,-3];

    let xs = [feetL,kneeL,feetR,kneeR,hip,torso,head,armL,armR];
    xs = micos.math_utils.translate(xs,[0,3]);

    //lets scale down the xs
    xs = micos.math_utils.scale(xs,0.7,0.7);

    let constraints_distance = [[0,1,0.7],[1,4,0.7],[2,3,0.7],[3,4,0.7],[4,5,0.7],[5,6,0.7],[5,7,0.7],[5,8,0.7]];

    let masses = [1,1,1,1,5,5,3,1,1];

    let dampingCoef = 0.07;

    //no ground

    let polygons = [];

    // a simple spring to point attached to the left foot and  to 0,0
    //springs2points must be like [[idx1,fixed_point, spring_constant, rest_length],...
    let springs2points = [[0,[0,4],5,1]];

    let STATE = {

        "xs": xs,
        "vs": xs.map((x) => {return [0,0]}),
        "masses": masses,
        "constraints_distance": constraints_distance,
        "constraints_pin": [],
        "polygons": polygons,
        "t":0,
        "gravity": [1,[0,-1]],
        "springs": [],//particle1, particle2, stiffness, restLength
        "damping": xs.map((x,i) => [[i],dampingCoef]),
        "springs2points": springs2points,

    };

    return STATE;

}


function spiderState(legLength){

    //a simple spider like state, a big mass in the head, and 2 legs comprised of 2 segments each
    // pin constraints at the feet so   the masses are feet->knee->head->knee->feet

    //at an initial angle of 60º

    let feetL = [0,0];
    let gradRad = 60*Math.PI/180;
    let kneeL = math.add(feetL,math.multiply(legLength, [Math.cos(gradRad), Math.sin(gradRad)]));
    let head = math.add(kneeL,math.multiply(legLength, [Math.cos(-gradRad), Math.sin(-gradRad)]));

    let kneeR = math.add(head, math.multiply(legLength, [Math.cos(gradRad), Math.sin(gradRad)]));
    let feetR = math.add(kneeR, math.multiply(legLength, [Math.cos(-gradRad), Math.sin(-gradRad)]));

    let xs = [kneeL,head,kneeR];
    //lets translate xs to the left
    xs = math_utils.translate(xs,[-0.7,0]);
    feetL = math.add(feetL,[-0.7,0]);
    feetR = math.add(feetR,[-0.7,0]);

    let masses = [1,8,1];
    let constraints_distance = [
                                [0,1,legLength],
                                [1,2,legLength],
                                ];

    let dampingCoef = 0.1;

    //let ground = [[-4,-1],[-4,-0.5],[4,-0.5],[4,-1]];

    let polygons = [];

    let STATE = {
        "xs": xs,
        "vs": xs.map((x) => {return [0,0]}),
        "masses": masses,
        "constraints_distance": constraints_distance,
        "constraints_pin": [[0,feetL,legLength],[2,feetR,legLength]],
        "polygons": polygons,
        "t":0,
        "gravity": [0.4,[0,-1]],
        "springs": [],//particle1, particle2, stiffness, restLength
        "damping": xs.map((x,i) => [[i],dampingCoef]),

    };

    return STATE;




}

function coupledDoublePendulum(L =1, g = 1, mass = 1){
    //two double pendulums coupled by a spring between the first masses

    let x00 = [0,-L];
    let x01 = [x00[0] , x00[1]-L];
    let x10 = [x01[0] + L, -L];
    let x11 = [x10[0], x10[1]-L];
    let xs = [x00,x01,x10,x11];
    //lets translate xs to the left

    let shift = -0.5;
    xs = math_utils.translate(xs,[shift,-shift]);

    let pinPoint0 = [xs[0][0],xs[0][1]+L];
    let pinPoint1 = [xs[2][0],xs[2][1]+L];


    let v0 = [[3,0],[0,0],[1,0],[0,0]];
    let masses = [mass,mass,mass,mass];
    let constraints_distance = [[0,1,L],[2,3,L]];
    let constraints_pin = [[0,pinPoint0,L],[2,pinPoint1,L]];
    let gravity = [g, [0,-1]];
    let time = 0;
    let dampingCoef = 0.05;
    let springs = [[0,2,4,L]];

    let initialState = {
        xs: xs,
        vs: v0,
        masses: masses,
        constraints_distance: constraints_distance,
        constraints_pin: constraints_pin,
        gravity: gravity,
        time: time,
        damping : xs.map((x,i) => [[i],dampingCoef]),
        springs: springs,
        }
    

    return initialState;


}





function rollingCircleSystemConfig(nparticles, radius,g, omega,mass = 1){
    //simple recreation of a rolling circle out of particles
    //to build the needed constraints so the circle is a rigid body we need to get neighs from delauney triangulation

    //lets create a xs with the n particles, we set the origin at (0,radius) and we place the particles in a circle in clockwise order
    
    //now for the vs we need to calculate the velocity of the center of mass and the angular velocity
    // rolling condition is v_cm = omega x r
    //velocity of a particle on the circle is v_cm+ omega x r_particle
    //rolling to the right, clockwise rotation means omega is negative so vparticle = omegaR \vec{e_{x}} -|omega|R\vec{\theta}
    //\vec{\theta} = (-sin(theta),cos(theta))  


    let xs = [];
    let vs = [];
    let origin = [0,radius];
    let angle = 2*Math.PI/nparticles;
    let vCom = [radius*omega,0];
    for (let i = 0; i < nparticles; i++){
        let anglei = i*angle
        let x = [radius*Math.cos(anglei), radius*Math.sin(anglei)];
        xs.push(x);
        let vang = [omega*x[1], -omega*x[0]];
        let vparticle = math.add(vCom,vang);
        vs.push(vparticle);
    }

    xs = micos.math_utils.translate(xs, origin);

    let constraints_distance = getConstraintsRigid(xs); //rigidification baby

    let masses = Array(nparticles).fill(mass);

    let polGround = [[-3,0],[3,0],[3,-1],[-3,-1]];

    let initialState = {
        xs: xs,
        vs: vs,
        masses: masses,
        constraints_distance: constraints_distance,
        constraints_pin: [],
        collisions: [],
        gravity: [g, [0,-1]],
        time: 0,
        polygons: [polGround],
        }

    return initialState;
    }




    function tumblingBoxSystemConfig(side, g, Lwedge, angle){
        //here we aim to simulate the tumbling box along a wedge
        let origin = [0,0];
        let xs = [[origin[0],origin[1]],
                  [origin[0],origin[1]+side],
                    [origin[0]+side,origin[1]+side],
                    [origin[0]+side,origin[1]]];
    
        let constraints_distance = getConstraintsRigid(xs)
        
        
    
        let base = Lwedge;
        //angle in deg
        let height = base*Math.tan(angle*Math.PI/180);
        let wedge = [[0,0],[0,height],[base,0]];
        //let wedge = [[0,0],[0,1],[1,0]];
        //now, lets place the box on the wedge, we'll place the pivot box vertex at intersection of the 90º angle of the wedge bisector and the wedge slope
        // for this we start by shifting the box by [-side,0] to the left so the pivot is at [0,0], we then rotate the box around the pivot by angle
        //and finally we shift the box by the vector slope_midpoint+origin. 
        //we'll use translate and rotate functions from math_utils.js
       // xs = translate(xs, [-side,0]);
        xs = rotate(xs,-angle,[0,0]);
        //midpoint [0,height]+[base/2,height]/2
        let slopeMidpoint = [base/2,height-height/2];
        wedge = micos.math_utils.translate(wedge, slopeMidpoint.map(x => -x));
    
        //xs = translate(xs, slopeMidpoint);
    
        let STATE = {
            xs: xs,
            constraints_distance: constraints_distance,
            constraints_pin: [],
            "constraints_values": [],
            "constraints_mae": 0,
            "targets": [],
            "masses": [1,1,1,1],
            "vs": [[0,0],[0,0],[0,0],[0,0]],
            "collisions": [],
            "J": null,
            "dotJ": null,
            "t":0,
            "gravity": [g,[0,-1]],
            "polygons": [wedge], 
            "springs2points": []
    
    
        };
    
        //lets translate box and wedge a little to the left and down
        return STATE;
    }
    


    function dumbellSlidingOnWallSystemConfig(L,mass1,mass2,angle, g){
        //simple test for contact constraints, without friction
        // a dumbell system, two masses with a distance constraint, mass1 is on the ground
        //mass2 is on the wall, angle is the initial angle of the dumbell with the horizontal
        //the system will slide down the wall with an effective 1 degree of freedom (the angle)
        //we must create two polygons, one for the ground and one for the wall
        // and we place mass1 at L(cos angle,0) and mass2 at L(0, sin angle)
    
        // we build polygons clock-wisely so the edges and normals are correctly computed
    
        let polWall = [[0,-3],[-1,-3],[-1,3],[0,3]];
        let polGround = [[-3,0],[3,0],[3,-1],[-3,-1]];
        angle = angle*Math.PI/180;
        let x1 = [L*Math.cos(angle),0];
        let x2 = [0,L*Math.sin(angle)];
        let x0 = [x1,x2];
        let v0 = [[0,0],[0,0]];
        let masses = [mass1,mass2];
    
        let initialState = {
            xs: x0,
            vs: v0,
            masses: masses,
            constraints_distance: [[0,1,L]],
            constraints_pin: [],
            collisions: [],
            gravity: [g,[0,-1]],
            time: 0,
            polygons: [polGround,polWall],
            }
    
    
    
        return initialState;
    }
    
    function dumbellSlidingOnWallAnalyticalSolution(L, mass, angle, g, tArray){
        // we'll integrate the angle evolution equation for the simplified assumption of equal masses
        //Lagrangian = T-V. T = 1/2 m1 v1^2 + 1/2 m2 v2^2   and V = mgx2 = mgLsin(theta)
        //v1 = d_t(Lcos(angle)) = -Lsin(angle)theta_dot v2 = d_t(Lsin(angle)) = Lcos(angle)theta_dot
        //T = 1/2 m1 L^2 sin^2(theta) theta_dot^2 + 1/2 m2 L^2 cos^2(theta) theta_dot^2 = mL^{2}theta_dot^2/2
        //V = mgLsin(theta)
        //L = mL^{2}theta_dot^2 - mgLsin(theta), eofm d_t(partial_L/partial_theta_dot) - partial_L/partial_theta = 0
        //2mL^{2}theta_dot theta_ddot + mgLcos(theta) = 0
        //theta_ddot = -g/(2L)cos(theta)
    
        let theta0 = angle*Math.PI/180;
    
        let dxdtFun = (x,t) => {
            let [theta,theta_dot] = x;
            let dthetadt = theta_dot;
            let dtheta_dotdt = (-g/(L))*Math.cos(theta);
            return [dthetadt,dtheta_dotdt];
        
        }
    
        let x0 = [theta0,0];
        let out = micos.math_utils.integrateOdeOnTarray(dxdtFun, x0, tArray,"Euler")//out [tarray.length,2] 
        let thetaArray = out.map(x => x[0]);
        let theta_dotArray = out.map(x => x[1]);
    
        let xsArray = thetaArray.map((theta,i) => [ [L*Math.cos(theta),0], [0,L*Math.sin(theta)] ]);
        let vsArray = theta_dotArray.map((theta_dot,i) => [ [-L*Math.sin(thetaArray[i])*theta_dot,0], [0,L*Math.cos(thetaArray[i])*theta_dot] ]);
        let externalForces = tArray.map(t => [ [0,-mass*g],[0,-mass*g] ]);
        let constraintForces = micos.solver.calculateConstraintForcesFromTrajectory(vsArray, [mass,mass], externalForces, tArray);
    
        let stateStory = [];
    
        for (let i = 0; i < tArray.length; i++){
            let state = {
                "xs": xsArray[i],
                "vs": vsArray[i],
                "time": tArray[i],
                "constraint_forces": constraintForces[i],
                "external_forces": externalForces[i],
                "gravity": [g,[0,-1]],
                "constraints_distance": [[0,1,L]],
                "polygons": [[[-3,0],[3,0],[3,-1],[-3,-1]],[[0,-3],[-1,-3],[-1,3],[0,3]]],
                "masses": [mass,mass],
    
            }
            stateStory.push(state);
        }
    
        return stateStory;
        
    }
    

    function pendulumSystemConfig(thetaO, omegaO, g, L, mass){

        let x0 = [[L*Math.sin(thetaO), -L*Math.cos(thetaO)]]; // [nparticles,2] [1,2] ( 1 particle)
        let v0 = [[omegaO*L*Math.cos(thetaO), omegaO*L*Math.sin(thetaO)]];
        let masses = [mass];
        let constraints_distance = [];
        let constraints_pin = [[0, [0,0], L]];
        let gravity = [g, [0,-1]];
        let time = 0;
        let initialState = {
            xs: x0,
            vs: v0,
            masses: masses,
            constraints_distance: constraints_distance,
            constraints_pin: constraints_pin,
            collisions: [],
            gravity: gravity,
            time: time,
            polygons: [],
            }
    
        
        
        return initialState;
    
    }

    function doublePendulumSystemConfig(thetaO1, omegaO1, thetaO2, omegaO2, g, L1, L2, mass1, mass2){
        let x00 = [L1*Math.sin(thetaO1), -L1*Math.cos(thetaO1)];
        let x01 = [x00[0] + L2*Math.sin(thetaO2), x00[1] - L2*Math.cos(thetaO2)];
        let x0 = [x00, x01];
        let v0 = [[omegaO1*L1*Math.cos(thetaO1), omegaO1*L1*Math.sin(thetaO1)], [omegaO2*L2*Math.cos(thetaO2), omegaO2*L2*Math.sin(thetaO2)]];
        let masses = [mass1, mass2];
        let constraints_distance = [ [0,1, L2]];
        let constraints_pin = [[0, [0,0], L1]];
        let gravity = [g, [0,-1]];
        let time = 0;
        let dampingCoef = 0.014;
        let initialState = {
            xs: x0,
            vs: v0,
            masses: masses,
            constraints_distance: constraints_distance,
            constraints_pin: constraints_pin,
            gravity: gravity,
            time: time,
            damping : x0.map((x,i) => [[i],dampingCoef]),
            }
    
        return initialState;
    }
    


    function chainSystemConfig(nparticles, g, L, mass, angle){
        // we'll create a chain of particles with a given angle to avoid an initial static configuration
        let x0 = [];
        let v0 = [];
        let masses = [];
        let constraints_distance = [];
        let constraints_pin = [];
        let gravity = [g, [0,-1]];
        let time = 0;
        constraints_pin.push([0, [0,0], L]);
        let dampingCoef = 0.03;
        for (let i = 0; i < nparticles; i++){
    
            // lets add a little angle perturbation to avoid singular configurations
            let angle_perturbation = angle;
            let x = [(i+1)*L*Math.cos(angle_perturbation), -(i+1)*L*Math.sin(angle_perturbation)];
            let v = [0,0];
            x0.push(x);
            v0.push(v);
            masses.push(mass);
            if (i > 0){
                constraints_distance.push([i-1,i,L]);
            }
        }
    
        let initialState = {
            xs: x0,
            vs: v0,
            masses: masses,
            constraints_distance: constraints_distance,
            constraints_pin: constraints_pin,
            gravity: gravity,
            time: time,
            damping : x0.map((x,i) => [[i],dampingCoef]),

            }
    
        return initialState;
    }



export {compassWalkerState, rollingCircleSystemConfig, tumblingBoxSystemConfig,
     pendulumSystemConfig, doublePendulumSystemConfig, chainSystemConfig, dumbellSlidingOnWallSystemConfig,
    spiderState, dumbellSlidingOnWallAnalyticalSolution, coupledDoublePendulum, bungeeJumperState, rigidDoublePendulum};
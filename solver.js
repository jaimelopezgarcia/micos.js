

function stepSemiEuler(STATE, PARAMETERS){
    /*
    *   Provisional semi-implicit euler integration
    */
    let xs = STATE.xs;
    let vs = STATE.vs;
    let masses = STATE.masses;
    let external_forces = STATE.external_forces;
    let constraint_forces = STATE.constraint_forces;
    let contact_forces = STATE.contact_forces;
    let friction_forces = STATE.friction_forces;
    let total_forces = xs.map((x,i) => math.add(external_forces[i], constraint_forces[i], contact_forces[i], friction_forces[i]));

    let dt = PARAMETERS["dt"];
    // gravity is assumed to be computed already in the external forces

    let new_vs = vs.map((v,i) => math.add(v, math.multiply(total_forces[i], dt/masses[i])));
    let new_xs = xs.map((x,i) => math.add(x, math.multiply(new_vs[i], dt)));
    let new_time = STATE.time + dt;


    let newState = {
        "xs": new_xs,
        "vs": new_vs,
        "time": new_time,
    }

    // lets add the rest of properties to newState copying them from STATE, excluding ["xs", "vs","time"]
    for (let key in STATE){
        if (!["xs", "vs", "time"].includes(key)){
            newState[key] = STATE[key];
        }
    }

    return newState;
}

export { stepSemiEuler };
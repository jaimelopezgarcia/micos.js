import {getSvgRelativeCoords} from "./utils.js";

const DEBUG = true;



function updateVisualOptions(key, value, VISUAL_OPTIONS){
    VISUAL_OPTIONS[key] = value;
}
function populateVisualOptionsPanel(VISUAL_OPTIONS, id_panel){
    /*
    *   Function to populate the visual options panel with the visual options
    *   @param VISUAL_OPTIONS: Object containing the visual options (differenciates between numerical and string options)
    *   @param id_panel: ID of the panel to populate
    *   @return: void
    */
    let visualOptionsPanel = document.getElementById(id_panel);
    visualOptionsPanel.appendChild(document.createElement("br"));
    for (let key in VISUAL_OPTIONS){
        let input = document.createElement("input");
        input.type = "text";
        input.id = id_panel + "-" + key;
        input.value = VISUAL_OPTIONS[key];
        input.onchange = function(){
            let val = input.value;
            if (!isNaN(parseFloat(val))){
                val = parseFloat(val);
            }
            updateVisualOptions(key, val, VISUAL_OPTIONS);
            if (DEBUG){
                console.log("Visual Option " + key + " updated to " + val);
            }
        }
        let label = document.createElement("label");
        label.htmlFor = input.id;
        label.innerHTML = key;
        
        let label_inputDiv = document.createElement("div");
        label_inputDiv.appendChild(label);
        label_inputDiv.appendChild(input);
        visualOptionsPanel.appendChild(label_inputDiv);
        label_inputDiv.classList.add("label-input");

    }
}

function updateParameters(key, value, PARAMETERS, OBSERVERS_PARAMETERS){
    PARAMETERS[key] = value;
    OBSERVERS_PARAMETERS.forEach(observer => observer(PARAMETERS));
}
function replaceParameters(new_parameters, PARAMETERS, OBSERVERS_PARAMETERS){
    for (let key in new_parameters){
        if (PARAMETERS[key] !== undefined){
            throw new Error("Parameter " + key + " is not defined in the PARAMETERS object");
        }
        PARAMETERS[key] = new_parameters[key];
    }
    OBSERVERS_PARAMETERS.forEach(observer => observer(PARAMETERS));
}



function populateParametersPanel(PARAMETERS, id_panel, OBSERVERS_PARAMETERS){
    /*
    *   Function to populate the parameters panel with the parameters
    *   @param PARAMETERS: Object containing the parameters
    *  @param id_panel: ID of the panel to populate
    *  @return: void
    * @note: The parameters are numerical
    */
    let parametersPanel = document.getElementById(id_panel);
    parametersPanel.appendChild(document.createElement("br"));
    

    for (let key in PARAMETERS){
        let input = document.createElement("input");
        input.type = "text";
        input.id = id_panel + "-" + key;
        input.value = PARAMETERS[key];
        input.onchange = function(){
            let val = parseFloat(input.value);
            updateParameters(key, val, PARAMETERS, OBSERVERS_PARAMETERS);
            if (DEBUG){
                console.log("Parameter " + key + " updated to " + val);
            }
        }
        let label = document.createElement("label");
        label.htmlFor = input.id;
        label.innerHTML = key;
        
        let label_inputDiv = document.createElement("div");
        label_inputDiv.appendChild(label);
        label_inputDiv.appendChild(input);
        parametersPanel.appendChild(label_inputDiv);
        label_inputDiv.classList.add("label-input");

    }
}






function makeParticleDraggable(SVG,particle_id, particle_idx,
                                STATE,updateState,canvas2Model, callbacks = [],
                                ){
    /*
    *   Function to make a particle draggable
    *  Assumes that the particle is a circle with a text element inside
    *  And the particle_id last part is the particle_idx
    * STATE must have a xs key with the positions of the particles
    * updateState is a function to update the state, intended to channel all STATE updates
    * canvas2Model must be a function to convert canvas to model coordinates
    * It must have signature canvas2Model(pointsCanvas) -> pointsModel arrays of [npoints, 2]
    */
    //lets throw an error if the number after the last _ in the id doesnt match the particle_idx
    let split_id = particle_id.split("_");
    let particle_number = split_id[split_id.length-1];
    if (particle_number != particle_idx){
        throw new Error(`particle_number ${particle_number} does not match particle_idx ${particle_idx}`);
    }
    let particle = document.getElementById(particle_id);
    //lets throw an error if the particle is not found
    if (particle === null){
        throw new Error("Particle not found with id " + particle_id);
    }
    particle.addEventListener("mousedown", function(event){
        event.preventDefault();
        let [xabs, yabs] = [event.clientX, event.clientY];
        let [xsvg, ysvg] = getSvgRelativeCoords(SVG, xabs, yabs);

        function mouseMoveCallback(event){
            let [xabs, yabs] = [event.clientX, event.clientY];
            let [xsvg, ysvg] = getSvgRelativeCoords(SVG, xabs, yabs);
            particle.querySelector("circle").setAttribute("cx", xsvg);
            particle.querySelector("circle").setAttribute("cy", ysvg);
            particle.querySelector("text").setAttribute("x", xsvg);
            particle.querySelector("text").setAttribute("y", ysvg);

            let [xnewCanvas, ynewCanvas] =[xsvg, ysvg]
            let [xnewModel, ynewModel] = canvas2Model([[xnewCanvas, ynewCanvas]])[0];
            let xs = STATE.xs;
            xs[particle_idx] = [xnewModel, ynewModel];
            updateState("xs", xs,STATE);
        }
        function mouseUpCallback(event){

            for (let callback of callbacks){
                callback();
            }
            
            if (DEBUG){
                console.log("UPDATED POSITIONS PARTICLE", particle_idx, JSON.stringify(STATE.xs));
            }   
            SVG.removeEventListener("mousemove", mouseMoveCallback);
            SVG.removeEventListener("mouseup", mouseUpCallback);
        }
        SVG.addEventListener("mousemove", mouseMoveCallback);
        SVG.addEventListener("mouseup", mouseUpCallback);
    })

}

function makeCircleDraggable(SVG, particle_id,
                                updateCallback, 
                                 canvas2Model,
                                  callbacks = []){
                /*
                *   A more general version of makeParticleDraggable
                * The updating logic is handled by the updateCallback function declared in the parent scope
                * It must have signature updateCallback(point) 
                */
                //lets throw an error if the number after the last _ in the id doesnt match the particle_idx
                let split_id = particle_id.split("_");
                let particle_idx = split_id[split_id.length-1];

                
                let particle = document.getElementById(particle_id);
                //lets throw an error if the particle is not found
                if (particle === null){
                throw new Error("Particle not found with id " + particle_id);
                }
                particle.addEventListener("mousedown", function(event){
                event.preventDefault();
                let [xabs, yabs] = [event.clientX, event.clientY];
                let [xsvg, ysvg] = getSvgRelativeCoords(SVG, xabs, yabs);

                function mouseMoveCallback(event){
                    let [xabs, yabs] = [event.clientX, event.clientY];
                    let [xsvg, ysvg] = getSvgRelativeCoords(SVG, xabs, yabs);
                    particle.querySelector("circle").setAttribute("cx", xsvg);
                    particle.querySelector("circle").setAttribute("cy", ysvg);
                    particle.querySelector("text").setAttribute("x", xsvg);
                    particle.querySelector("text").setAttribute("y", ysvg);

                    let [xnewCanvas, ynewCanvas] =[xsvg, ysvg]
                    let [xnewModel, ynewModel] = canvas2Model([[xnewCanvas, ynewCanvas]])[0];
                        
                    let point = [xnewModel, ynewModel];
                    updateCallback(point);
                }


                function mouseUpCallback(event){

                for (let callback of callbacks){
                callback();
                }

      
                SVG.removeEventListener("mousemove", mouseMoveCallback);
                SVG.removeEventListener("mouseup", mouseUpCallback);
                }
                SVG.addEventListener("mousemove", mouseMoveCallback);
                SVG.addEventListener("mouseup", mouseUpCallback);
                })

}

export {populateVisualOptionsPanel, populateParametersPanel,
        updateParameters, replaceParameters,updateVisualOptions,makeParticleDraggable, makeCircleDraggable};
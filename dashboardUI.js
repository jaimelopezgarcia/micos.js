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










export {populateVisualOptionsPanel, populateParametersPanel,
        updateParameters, replaceParameters,updateVisualOptions};
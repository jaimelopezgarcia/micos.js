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


function particleInteractionMouse(SVG, particle_id, canvas2Model, mouseMoveCallbacks = [], mouseUpCallbacks = []){
    //this function is intended to execute the provided callbacks when the particle is clicked
    // all mouseMoveCalbacks must have the signature callback(mouseCoordsModel, particle_idx)
    // all mouseUpCallbacks must have the signature callback()
    // with mouseup we remove the mousemove and mouseup event listeners
    let particle = document.getElementById(particle_id);
    //lets throw an error if the particle is not found
    if (particle === null){
        throw new Error("Particle not found with id " + particle_id);
    }
    //lets retrieve the particle_idx ( it must be the last string after _ in the id)
    let split_id = particle_id.split("_");
    //parseInt
    let particle_idx = parseInt(split_id[split_id.length-1]);
    //if not a number, throw an error
    if (isNaN(particle_idx)){
        throw new Error("Particle id must end with a number");
    }
    particle.addEventListener("mousedown", function(event){
        event.preventDefault();

        function mouseMoveCallback(event){
            let [xabs, yabs] = [event.clientX, event.clientY];
            let [xsvg, ysvg] = getSvgRelativeCoords(SVG, xabs, yabs);
            let [xnewCanvas, ynewCanvas] =[xsvg, ysvg]
            let [xnewModel, ynewModel] = canvas2Model([[xnewCanvas, ynewCanvas]])[0];
            let mouseCoordsModel = [xnewModel, ynewModel];
            for (let callback of mouseMoveCallbacks){
                callback(mouseCoordsModel, particle_idx);
            }
        }
        function mouseUpCallback(event){
            SVG.removeEventListener("mousemove", mouseMoveCallback);
            SVG.removeEventListener("mouseup", mouseUpCallback);

            for (let callback of mouseUpCallbacks){
                callback();
            }
        }

        SVG.addEventListener("mousemove", mouseMoveCallback);
        SVG.addEventListener("mouseup", mouseUpCallback);

    })

}
       





function makeParticleDraggable(SVG,particle_id,
                                STATE,updateState,
                                canvas2Model, callbacks = [],
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
    let particle_idx = parseInt(particle_number);
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


class SelectionPanel {
    /*
    Example usage
    const data = {
        Particle_0: ['X Position', 'Y Position'],
        Particle_1: ['X Position', 'Y Position']
    };

    const selectionPanel = new SelectionPanel(data);

    //methods
    selectionPanel.gatherSelections(); //returns an array of the selected checkboxes
    selectionPanel.clearSelections(); //clears all the selected checkboxes
    selectionPanel.selectAll(); //selects all the checkboxes
    

    Example of the html structure that will be created
        <div class="selectionPanel">
            <ul>
                <li>
                    <span>Particle_0</span>
                    <ul>
                        <li><input type="checkbox" id="particle_0_x"> X Position</li>
                        <li><input type="checkbox" id="particle_0_y"> Y Position</li>
                    </ul>
                </li>
                <li>
                    <span>Particle_1</span>
                    <ul>
                        <li><input type="checkbox" id="particle_1_x"> X Position</li>
                        <li><input type="checkbox" id="particle_1_y"> Y Position</li>
                    </ul>
                </li>
            </ul>
    */
    constructor(data, parentDivId, id = "selectionPanel") {
        this.id = id;
        const selector = `#${id}`;

        this.selector = selector;
        this.data = data;
        this.parentDiv = document.getElementById(parentDivId);
        if (!this.parentDiv) {
            console.error('Parent div not found');
            return;
        }
        this.selectionPanelDiv = this.init();
    }

    init() {
        let selectionPanelDiv = this.createSelectionPanel(this.data, this.parentDiv);
        this.injectStyles();
        this.attachEventListeners();

        return selectionPanelDiv;
    }

    getPanel() {
        return this.selectionPanelDiv;
    }

    createSelectionPanel(data,parentDiv) {

        let selectionPanel = document.querySelector(`${this.selector}`);
        // if not null we throw an error of non-unique id
        if (selectionPanel) {
            console.error('Non-unique id');
            return;
        }

        selectionPanel = document.createElement('div');
        selectionPanel.id = this.id;
        parentDiv.appendChild(selectionPanel);
    

        //lets assert the correct format of data structure  {key: [valueString1, valueString2, ...]}
        if (!data || typeof data !== 'object') {
            console.error('Invalid data format');
            return;
        }


            
            let ul = document.createElement('ul');
            for (let key in data) {
                let li = document.createElement('li');
                let span = document.createElement('span');
                span.textContent = key;
                li.appendChild(span);
                let subUl = document.createElement('ul');
                for (let i = 0; i < data[key].length; i++) {
                    let subLi = document.createElement('li');
                    let input = document.createElement('input');
                    input.type = 'checkbox';
                    //for id we use key_valueStringIndex
                    input.id = `${key}->${data[key][i]}`;
                    let label = document.createElement('label');
                    label.htmlFor = input.id;
                    label.textContent = data[key][i];
                    subLi.appendChild(input);
                    subLi.appendChild(label);
                    subUl.appendChild(subLi);
                }
                li.appendChild(subUl);
                ul.appendChild(li);
            }

            selectionPanel.appendChild(ul);

        return selectionPanel;
        }

    injectStyles() {
        const styles = `
            ${this.selector} ul {
                list-style-type: none;
            }

            ${this.selector} li > ul {
                display: none;
                margin-left: 20px;
            }
            ${this.selector} li > span:hover + ul,
            ${this.selector} li > span:focus + ul {
                display: block;
            }
            ${this.selector} input[type="checkbox"]:checked + label {
                font-weight: bold;
            }

            ${this.selector} li > span::before {
            content: '▶';
            display: inline-block;
            margin-right: 5px;
            transform: rotate(0deg);
            transition: transform 0.3s ease;
        }

        /* Rotate marker when ul is displayed */
        ${this.selector} li > span:hover::before,
        ${this.selector} li > span:focus::before {
            transform: rotate(90deg);
        }
        `;
        //lets append it to the existing style tag
        let styleTag = document.querySelector('style');
        if (!styleTag) {
            styleTag = document.createElement('style');
            document.head.appendChild(styleTag);
        }
        styleTag.textContent += styles;

    }

    attachEventListeners() {
        document.querySelectorAll(`${this.selector} li > span`).forEach(item => {
            item.addEventListener('click', () => {
                let subList = item.nextElementSibling;
                subList.style.display = subList.style.display === "block" ? "none" : "block";
            });
        });
    }

    gatherSelections() {
        let selections = [];
        document.querySelectorAll(`${this.selector} input[type="checkbox"]:checked`).forEach(item => {
            selections.push(item.id);
        });
        return selections;
    }

    clearSelections() {
        document.querySelectorAll(`${this.selector} input[type="checkbox"]:checked`).forEach(item => {
            item.checked = false;
        });
    }

    selectAll() {
        document.querySelectorAll(`${this.selector} input[type="checkbox"]`).forEach(item => {
            item.checked = true;
        });
    }


}


export {populateVisualOptionsPanel, populateParametersPanel,
        updateParameters, replaceParameters,updateVisualOptions,makeParticleDraggable, makeCircleDraggable, particleInteractionMouse,
        SelectionPanel
        };
import {Drawer} from "./drawing.js";
import {step} from "./solver.js";


function play(svg,STATE,
    PARAMETERS = {},callbacksActuators = [], msBetweenFrames = 5){

    let drawer = new Drawer(svg);
    drawer.drawState(STATE);               
    let mainCallback = () => {
    step(STATE, PARAMETERS)
    drawer.drawState(STATE);
    }


    if (window.intervalId!== null){
    clearInterval(window.intervalId);
    }

    window.intervalId = setInterval(mainCallback, msBetweenFrames);


    let isPaused = false;
    document.addEventListener("keydown", (event) => {
    if (event.code === "Space"){
    if (isPaused){
        window.intervalId = setInterval(mainCallback, msBetweenFrames);
    }else{
        clearInterval(window.intervalId);
    }
    isPaused = !isPaused;  
    }
});


}



export {play};
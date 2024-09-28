import {Drawer} from "./drawing.js";
import {step} from "./solver.js";

class Profiler {
    constructor(smoothingFactor = 0.1) {
        this.timers = {}; // Holds start times for each named timer
        this.records = {}; // Holds stats for each named timer (last_time, avg, max_time)
        this.smoothingFactor = smoothingFactor; // Smoothing factor for exponential average
    }

    // Start or reset a timer for the given name
    _create_timer(name) {

        //it only initializes the first time

        this.timers[name] = performance.now(); // Set current time for the given name
        // Initialize records if they don't exist for this timer
        if (!this.records[name]) {
            this.records[name] = { last_time: 0, avg: 0, max_time: 0 };
        }
    }

    start_timer(name) {
        if (!this.timers[name]) {
            console.info(`Timer '${name}' not found. Creating new timer.`);
            this._create_timer(name);
        }
        this.timers[name] = performance.now();

    }
    // Update and record the time elapsed since last set_timer or time() call
    time(name) {
        if (!this.timers[name]) {
            console.warn(`Timer '${name}' is not set.`);
            return;
        }

        const currentTime = performance.now();
        const elapsedTime = currentTime - this.timers[name];
        
        // Update timer for next call
        this.timers[name] = currentTime;

        // Update last_time
        this.records[name].last_time = elapsedTime;

        // Update smoothed running average using exponential smoothing
        this.records[name].avg = this.smoothingFactor * elapsedTime + (1 - this.smoothingFactor) * this.records[name].avg;

        // Update max_time
        if (elapsedTime > this.records[name].max_time) {
            this.records[name].max_time = elapsedTime;
        }
    }

    // Retrieve the timing statistics for the given name
    get_times(name) {
        if (!this.records[name]) {
            console.warn(`No records found for timer '${name}'.`);
            return null;
        }

        return {
            last_time: this.records[name].last_time,
            avg: this.records[name].avg,
            max_time: this.records[name].max_time
        };
    }
}



class Player {
    constructor(svg) {
        this.svg = svg;
        this.isPaused = false;
        this.isInit = false;
        this.opts = {};
        this.STATE = null;
        this.drawer = null;
        this.animationId = null;  // Store the animation ID to cancel the loop
    }

    _playInit() {
        this.isInit = true;

        if (this.isPaused) {
            return;
        }

        let lastTime = performance.now();
        let simTimeAccumulator = 0;
        let renderAccumulator = 0;  // Track time for rendering

        const gameLoop = () => {
            // Check if paused
            if (this.isPaused) return;

            // Access dynamic values directly from the `Player` instance (this.STATE, this.opts)
            const dt = this.opts.dt;
            const frameDuration = 1 / this.opts.targetFPS;
            const ratioSimTimeToRealTime = this.opts.ratioSimTimeToRealTime;

            const currentTime = performance.now();
            const realTimeDelta = (currentTime - lastTime) / 1000;  // Convert to seconds
            lastTime = currentTime;

            simTimeAccumulator += realTimeDelta * ratioSimTimeToRealTime;
            renderAccumulator += realTimeDelta;  // Accumulate time for rendering

            // Run the simulation steps as needed
            while (simTimeAccumulator >= dt) {
                step(this.STATE, { "dt": dt });  // Dynamic STATE reference
                simTimeAccumulator -= dt;
            }

            // Only render when enough time has passed
            if (renderAccumulator >= frameDuration) {
                this.drawer.drawState(this.STATE);  // Dynamic STATE reference for rendering
                renderAccumulator -= frameDuration;  // Subtract frameDuration instead of resetting
            }

            // Schedule the next loop iteration
            this.animationId = requestAnimationFrame(gameLoop);
        };

        // Start the game loop by calling it for the first time
        this.animationId = requestAnimationFrame(gameLoop);
    }

    play(STATE, options = {}) {
        this.STATE = STATE;

        let defaultOptions = {
            "dt": 0.005,
            "ratioSimTimeToRealTime": 2,
            "targetFPS": 60,
            "debug": false,
        };

        // Merge user-provided options with default options
        for (let key in options) {
            if (defaultOptions.hasOwnProperty(key)) {
                defaultOptions[key] = options[key];
            } else {
                console.error(`play: invalid option ${key}, valid options are ${Object.keys(defaultOptions)}`);
            }
        }

        this.opts = defaultOptions;

        if (this.drawer != null) {
            this.drawer.remove();  // Clear the previous drawer
            while (this.svg.firstChild) {  // Clear the SVG element
                this.svg.removeChild(this.svg.firstChild)
                }
        }

        this.drawer = new Drawer(this.svg, { "debug": this.opts.debug });

        // Only run _playInit once, loop will dynamically access updated attributes
        if (!this.isInit) {
            this._playInit();
        }
    }

    pausePlay() {
        this.isPaused = !this.isPaused;

        // Resume the loop if it was paused
        if (!this.isPaused && !this.animationId) {
            this._playInit();  // Restart the loop if unpaused
        }
    }
}

function plotstateStory(stateStory){

    let stateUtils = new StateUtils();//we'll use the getNormalizedConstraintViolations(STATE) to plot the const errors

    let layout = {
        title: 'State Story',
        xaxis: {
          title: 'Time'
        },
        yaxis: {
          title: 'Value'
        }
      };

    let i = 0;
    let xi = stateStory.map(state => state.xs[i][0]);
    let yi = stateStory.map(state =>  state.xs[i][1]);
    let vxi = stateStory.map(state =>  state.vs[i][0]);
    let vyi = stateStory.map(state =>  state.vs[i][1]);


    let time = stateStory.map(state =>  state.time);

    let vars = [xi,yi,vxi,vyi];
    let names = ["x","y","vx","vy"];

    let traces = [];
    for (let i = 0; i<vars.length;i++){

        let xdata = time;
        let ydata = vars[i];

        traces.push({
            x: xdata,
            y: ydata,
            mode: 'lines',
            name: names[i]
        });
        
    }


    let plotDiv = document.getElementById("plotlyDiv1");
    if (plotDiv==null){
        plotDiv = document.createElement("div");
        plotDiv.id = "plotlyDiv1";
        document.body.appendChild(plotDiv);
    }

    Plotly.newPlot('plotlyDiv1', traces, layout);



    let constErrors = stateStory.map(state => stateUtils.getNormalizedConstraintViolations(state));
    //constErrors is a nsteps x nconstraints array we have to create nconstraints traces
    let nConst = constErrors[0].length;
    let tracesConst = [];

    for (let i = 0;i<nConst;i++){
        let trace = constErrors.map((errors) => errors[i]);
        let xdata = time;
        let ydata = trace;

        tracesConst.push({
            x: xdata,
            y: ydata,
            mode: 'lines',
            name: "Constraint "+i
        });



        }

    let plotDivConst = document.getElementById("plotlyDiv2");
    if (plotDivConst==null){
        plotDivConst = document.createElement("div");
        plotDivConst.id = "plotlyDiv2";
        document.body.appendChild(plotDivConst);
    }

    Plotly.newPlot('plotlyDiv2', tracesConst, layout);

    }



export {Player};
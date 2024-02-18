const registerServiceWorker = async () => {
    if ("serviceWorker" in navigator) {
      try {
        const registration = await navigator.serviceWorker.register("/sw.js", {
          scope: "/",
        });
        if (registration.installing) {
          console.log("Service worker installing");
        } else if (registration.waiting) {
          console.log("Service worker installed");
        } else if (registration.active) {
          console.log("Service worker active");
        }
      } catch (error) {
        console.error(`Registration failed with ${error}`);
      }
    }
  };

registerServiceWorker();


function calculate_gravitational_force() {
    var mass1 = parseFloat(document.getElementById("mass1").value);
    var mass2 = parseFloat(document.getElementById("mass2").value);
    var distance = parseFloat(document.getElementById("distance").value);
  
    if (isNaN(mass1) || isNaN(mass2) || isNaN(distance)) {
      alert("Please enter valid numbers for mass and distance.");
      return;
    }
  
    var result = (6.67e-11 * mass1 * mass2) / Math.pow(distance, 2);
  
    var resultString = result.toString();
  
    document.getElementById("result").textContent = "Result: " + resultString;
  }

  function calculateIntensityRatio() {
    var mA = parseFloat(document.getElementById("mA").value);
    var mB = parseFloat(document.getElementById("mB").value);
    var iA = parseFloat(document.getElementById("iA").value);
    var iB = parseFloat(document.getElementById("iB").value);

    if (isNaN(mA) || isNaN(mB) || isNaN(iA) || isNaN(iB)) {
        alert("Please enter valid numbers for mA, mB, iA, and iB.");
        return;
    }

    var result;
    if (mA === Math.E) {
        result = -1 * ((2.5 * Math.log(iA / iB)) - mB);
    } else if (mB === Math.E) {
        result = (2.5 * Math.log(iA / iB)) + mA;
    } else if (iA === Math.E) {
        result = (10 * Math.pow((mA - mB), 0.4)) / iB;
    } else if (iB === Math.E) {
        result = iA / (10 * Math.pow((mA - mB), 0.4));
    }

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateSmallAngle() {
    var D = parseFloat(document.getElementById("D").value);
    var d = parseFloat(document.getElementById("d").value);

    if (isNaN(D) || isNaN(d)) {
        alert("Please enter valid numbers for D and d.");
        return;
    }

    var result = 206265 * (D / d);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateCircularVelocity() {
    var M = parseFloat(document.getElementById("M").value);
    var r = parseFloat(document.getElementById("r").value);

    if (isNaN(M) || isNaN(r)) {
        alert("Please enter valid numbers for M and r.");
        return;
    }

    var G = 6.67e-11;
    var result = Math.sqrt(G * M / r);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateResolvingPower() {
    var lam = parseFloat(document.getElementById("lam").value);
    var d = parseFloat(document.getElementById("d").value);

    if (isNaN(lam) || isNaN(d)) {
        alert("Please enter valid numbers for λ and d.");
        return;
    }

    var result = 251643 * lam / d;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateCompareLGP() {
    var DA = parseFloat(document.getElementById("DA").value);
    var DB = parseFloat(document.getElementById("DB").value);

    if (isNaN(DA) || isNaN(DB)) {
        alert("Please enter valid numbers for DA and DB.");
        return;
    }

    var result = Math.pow(DA / DB, 2);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateMagnification() {
    var FO = parseFloat(document.getElementById("FO").value);
    var FE = parseFloat(document.getElementById("FE").value);

    if (isNaN(FO) || isNaN(FE)) {
        alert("Please enter valid numbers for FO and FE.");
        return;
    }

    var result = FO / FE;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateWienLaw() {
    var T = parseFloat(document.getElementById("T").value);

    if (isNaN(T)) {
        alert("Please enter a valid number for T.");
        return;
    }

    var result = 2900000 / T;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateStefanBoltzmannLaw() {
    var T = parseFloat(document.getElementById("T").value);

    if (isNaN(T)) {
        alert("Please enter a valid number for T.");
        return;
    }

    var sigma = 5.67e-8;
    var result = sigma * Math.pow(T, 4);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateStefanBoltzmannLawLuminosity() {
    var R = parseFloat(document.getElementById("R").value);
    var T = parseFloat(document.getElementById("T").value);

    if (isNaN(R) || isNaN(T)) {
        alert("Please enter valid numbers for R and T.");
        return;
    }

    var sigma = 5.67e-8;
    var result = 4 * Math.PI * Math.pow(R, 2) * Math.pow(T, 4) * sigma;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateRedshift() {
    var lambdaO = parseFloat(document.getElementById("lambdaO").value);
    var lambdaE = parseFloat(document.getElementById("lambdaE").value);

    if (isNaN(lambdaO) || isNaN(lambdaE)) {
        alert("Please enter valid numbers for λO and λE.");
        return;
    }

    var result = lambdaO / lambdaE - 1;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateDopplerFormula() {
    var z = parseFloat(document.getElementById("z").value);

    if (isNaN(z)) {
        alert("Please enter a valid number for z.");
        return;
    }

    var c = 3e8;
    var result = c * z;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateRelativisticDopplerShift() {
    var v = parseFloat(document.getElementById("v").value);

    if (isNaN(v)) {
        alert("Please enter a valid number for v.");
        return;
    }

    var c = 3e8;
    var result = Math.sqrt((1 + v / c) / (1 - v / c));

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateHubbleLaw() {
    var d = parseFloat(document.getElementById("d").value);

    if (isNaN(d)) {
        alert("Please enter a valid number for d.");
        return;
    }

    var H0 = 70;
    var result = H0 * d;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateFusionEnergy() {
    var m = parseFloat(document.getElementById("m").value);

    if (isNaN(m)) {
        alert("Please enter a valid number for m.");
        return;
    }

    var c = 3e8;
    var result = m * Math.pow(c, 2);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateParallax() {
    var p = parseFloat(document.getElementById("p").value);

    if (isNaN(p)) {
        alert("Please enter a valid number for p.");
        return;
    }

    var result = 1 / p;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateFRatio() {
    var Lf = parseFloat(document.getElementById("Lf").value);
    var DO = parseFloat(document.getElementById("DO").value);

    if (isNaN(Lf) || isNaN(DO)) {
        alert("Please enter valid numbers for Lf and DO.");
        return;
    }

    var result = Lf / DO;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateDistanceModulus() {
    var m = parseFloat(document.getElementById("m").value);
    var M = parseFloat(document.getElementById("M").value);
    var d = parseFloat(document.getElementById("d").value);

    if (isNaN(m) || isNaN(M) || isNaN(d)) {
        alert("Please enter valid numbers for m, M, and d.");
        return;
    }

    var result = m - M + 5 - 5 * Math.log10(d);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateOrbitEccentricity() {
    var c = parseFloat(document.getElementById("c").value);
    var a = parseFloat(document.getElementById("a").value);

    if (isNaN(c) || isNaN(a)) {
        alert("Please enter valid numbers for c and a.");
        return;
    }

    var result = c / a;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateOrbitAphelion() {
    var a = parseFloat(document.getElementById("a").value);
    var e = parseFloat(document.getElementById("e").value);

    if (isNaN(a) || isNaN(e)) {
        alert("Please enter valid numbers for a and e.");
        return;
    }

    var result = a * (1 + e);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateOrbitPerihelion() {
    var a = parseFloat(document.getElementById("a").value);
    var e = parseFloat(document.getElementById("e").value);

    if (isNaN(a) || isNaN(e)) {
        alert("Please enter valid numbers for a and e.");
        return;
    }

    var result = a * (1 - e);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateKeplerLaw() {
    var a = parseFloat(document.getElementById("a").value);
    var p = parseFloat(document.getElementById("p").value);

    if (isNaN(a) || isNaN(p)) {
        alert("Please enter valid numbers for a and p.");
        return;
    }

    var result = Math.pow(a, 3) / Math.pow(p, 2);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateMassLuminosity() {
    var M = parseFloat(document.getElementById("M").value);

    if (isNaN(M)) {
        alert("Please enter a valid number for M.");
        return;
    }

    var result = Math.pow(M, 3.5);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateStellarLifetime() {
    var M = parseFloat(document.getElementById("M").value);

    if (isNaN(M)) {
        alert("Please enter a valid number for M.");
        return;
    }

    var result = 1e10 * Math.pow(M, -2.5);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateSchwarzschildRadius() {
    var m = parseFloat(document.getElementById("m").value);

    if (isNaN(m)) {
        alert("Please enter a valid number for m.");
        return;
    }

    var G = 6.67e-11;
    var c = 3e8;
    var result = (2 * G * m) / Math.pow(c, 2);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateUniversalGravitation() {
    var mA = parseFloat(document.getElementById("mA").value);
    var mB = parseFloat(document.getElementById("mB").value);
    var r = parseFloat(document.getElementById("r").value);

    if (isNaN(mA) || isNaN(mB) || isNaN(r)) {
        alert("Please enter valid numbers for mA, mB, and r.");
        return;
    }

    var G = 6.67e-11;
    var result = (G * mA * mB) / Math.pow(r, 2);

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateInverseSquares() {
    var L = parseFloat(document.getElementById("L").value);
    var d = parseFloat(document.getElementById("d").value);

    if (isNaN(L) || isNaN(d)) {
        alert("Please enter valid numbers for L and d.");
        return;
    }

    var result = L / (4 * Math.PI * Math.pow(d, 2));

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}

function calculateVacuumFrequency() {
    var lambda_val = parseFloat(document.getElementById("lambda_val").value);

    if (isNaN(lambda_val)) {
        alert("Please enter a valid number for lambda_val.");
        return;
    }

    var c = 3e8;
    var result = c / lambda_val;

    var resultString = result.toString();

    document.getElementById("result").textContent = "Result: " + resultString;
}




function openModal(contentId) {
    var modal = document.getElementById("myModal");
  var contentPlaceholder = document.getElementById("modal-content-placeholder");

  contentPlaceholder.innerHTML = "";
  switch (contentId) {
    case 'calculate_gravitational_force':
      contentPlaceholder.innerHTML = `
      <div class= "modal-container-user">
        <h2>Calculate Gravitational Force
        
        </h2>
        <div>
          <label for="mass1">Mass 1 (kg):</label>
          <input type="number" id="mass1">
        </div>
        <div>
          <label for="mass2">Mass 2 (kg):</label>
          <input type="number" id="mass2">
        </div>
        <div>
          <label for="distance">Distance (km):</label>
          <input type="number" id="distance">
        </div>
        <div>
          <label id ="result" class = "result">Result :</label>
        </div>
        <button onclick="calculate_gravitational_force()">Calculate</button>
        </div>
      `;
      break;
    case 'intensity_ratio':
      contentPlaceholder.innerHTML = `
        <h2>Intensity Ratio</h2>
        <div>
          <label for="mA">mA:</label>
          <input type="number" id="mA">
        </div>
        <div>
          <label for="mB">mB:</label>
          <input type="number" id="mB">
        </div>
        <div>
          <label for="dA">dA:</label>
          <input type="number" id="dA">
        </div>
        <div>
          <label for="dB">dB:</label>
          <input type="number" id="dB">
        </div>
        <div>
          <label id="result">Result:</label>
        </div>

        <button onclick="calculateIntensityRatio()">Calculate</button>
      `;
      break;
      case 'small_angle':
        contentPlaceholder.innerHTML = `
        <div class="modal-container-user">
        <h2>Calculate Small Angle</h2>
        <div>
          <label for="D">D:</label>
          <input type="number" id="D">
        </div>
        <div>
          <label for="d">d:</label>
          <input type="number" id="d">
        </div>
        <div>
          <label id="result" class="result">Result:</label>
        </div>
        <button onclick="calculateSmallAngle()">Calculate</button>
      </div>
        `;
        break;
        case 'circular_velocity':
            contentPlaceholder.innerHTML = `
            <div class="modal-container-user">
                <h2>Calculate Circular Velocity</h2>
                <div>
                    <label for="M">M:</label>
                    <input type="number" id="M">
                </div>
                <div>
                    <label for="r">r:</label>
                    <input type="number" id="r">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateCircularVelocity()">Calculate</button>
                </div>
            `;
            break;
            case 'resolving_power':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Resolving Power</h2>
                <div>
                    <label for="lam">λ (wavelength):</label>
                    <input type="number" id="lam">
                </div>
                <div>
                    <label for="d">d (diameter):</label>
                    <input type="number" id="d">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateResolvingPower()">Calculate</button>
                </div>
                    `;
                    break;
            case 'compare_LGP':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Compare LGP</h2>
                <div>
                    <label for="DA">DA:</label>
                    <input type="number" id="DA">
                </div>
                <div>
                    <label for="DB">DB:</label>
                    <input type="number" id="DB">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateCompareLGP()">Calculate</button>
                </div>
                `;
                break;
            case 'magnification':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Magnification</h2>
                <div>
                    <label for="FO">FO:</label>
                    <input type="number" id="FO">
                </div>
                <div>
                    <label for="FE">FE:</label>
                    <input type="number" id="FE">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateMagnification()">Calculate</button>
                </div>
                `;
                break;
            case 'wien_law':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Wien's Law</h2>
                <div>
                    <label for="T">T:</label>
                    <input type="number" id="T">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateWienLaw()">Calculate</button>
                </div>
                `;
                break;
            case 'stefan_boltzmann_law':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Stefan-Boltzmann Law</h2>
                <div>
                    <label for="T">T:</label>
                    <input type="number" id="T">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateStefanBoltzmannLaw()">Calculate</button>
                </div>
                `;
                break;
            case 'stefan_boltzmann_law_luminosity':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Stefan-Boltzmann Law Luminosity</h2>
                <div>
                    <label for="R">R:</label>
                    <input type="number" id="R">
                </div>
                <div>
                    <label for="T">T:</label>
                    <input type="number" id="T">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateStefanBoltzmannLawLuminosity()">Calculate</button>
                </div>
                `;
                break;
            case 'redshift':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Redshift</h2>
                <div>
                    <label for="lambdaO">λO:</label>
                    <input type="number" id="lambdaO">
                </div>
                <div>
                    <label for="lambdaE">λE:</label>
                    <input type="number" id="lambdaE">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateRedshift()">Calculate</button>
                </div>
                `;
                break;
            case 'doppler_formula':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Doppler Formula</h2>
                <div>
                    <label for="z">z:</label>
                    <input type="number" id="z">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateDopplerFormula()">Calculate</button>
                </div>
                `;
                break;
            case 'relativistic_doppler_shift':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Relativistic Doppler Shift</h2>
                <div>
                    <label for="v">v:</label>
                    <input type="number" id="v">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateRelativisticDopplerShift()">Calculate</button>
                </div>
                `;
                break;
            case 'hubble_law':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Hubble Law</h2>
                <div>
                    <label for="d">d:</label>
                    <input type="number" id="d">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateHubbleLaw()">Calculate</button>
                </div>
                `;
                break;
            case 'fusion_energy':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Fusion Energy</h2>
                <div>
                    <label for="m">m:</label>
                    <input type="number" id="m">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateFusionEnergy()">Calculate</button>
                </div>
                `;
                break;
            case 'parallax':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Parallax</h2>
                <div>
                    <label for="p">p:</label>
                    <input type="number" id="p">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateParallax()">Calculate</button>
                </div>
                `;
                break;
            case 'f_ratio':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate F Ratio</h2>
                <div>
                    <label for="Lf">Lf:</label>
                    <input type="number" id="Lf">
                </div>
                <div>
                    <label for="DO">DO:</label>
                    <input type="number" id="DO">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateFRatio()">Calculate</button>
                </div>
                `;
                break;
            case 'distance_modulus':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Distance Modulus</h2>
                <div>
                    <label for="m">m:</label>
                    <input type="number" id="m">
                </div>
                <div>
                    <label for="M">M:</label>
                    <input type="number" id="M">
                </div>
                <div>
                    <label for="d">d:</label>
                    <input type="number" id="d">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateDistanceModulus()">Calculate</button>
                </div>
                `;
                break;
            case 'orbit_eccentricity':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Orbit Eccentricity</h2>
                <div>
                    <label for="c">c:</label>
                    <input type="number" id="c">
                </div>
                <div>
                    <label for="a">a:</label>
                    <input type="number" id="a">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateOrbitEccentricity()">Calculate</button>
                </div>
                `;
                break;
            case 'orbit_aphelion':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Orbit Aphelion</h2>
                <div>
                    <label for="a">a:</label>
                    <input type="number" id="a">
                </div>
                <div>
                    <label for="e">e:</label>
                    <input type="number" id="e">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateOrbitAphelion()">Calculate</button>
                </div>
                `;
                break;
            case 'orbit_perihelion':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Orbit Perihelion</h2>
                <div>
                    <label for="a">a:</label>
                    <input type="number" id="a">
                </div>
                <div>
                    <label for="e">e:</label>
                    <input type="number" id="e">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateOrbitPerihelion()">Calculate</button>
                </div>
                `;
                break;
            case 'kepler_law':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Kepler's Law</h2>
                <div>
                    <label for="a">a:</label>
                    <input type="number" id="a">
                </div>
                <div>
                    <label for="p">p:</label>
                    <input type="number" id="p">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateKeplerLaw()">Calculate</button>
                </div>
                `;
                break;
            case 'mass_luminosity':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Mass-Luminosity</h2>
                <div>
                    <label for="M">M:</label>
                    <input type="number" id="M">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateMassLuminosity()">Calculate</button>
                </div>
                `;
                break;
            case 'stellar_lifetime':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Stellar Lifetime</h2>
                <div>
                    <label for="M">M:</label>
                    <input type="number" id="M">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateStellarLifetime()">Calculate</button>
                </div>
                `;
                break;
            case 'schwarzschild_radius':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Schwarzschild Radius</h2>
                <div>
                    <label for="m">m:</label>
                    <input type="number" id="m">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateSchwarzschildRadius()">Calculate</button>
                </div>
                `;
                break;
            case 'universal_gravitation':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Universal Gravitation</h2>
                <div>
                    <label for="mA">mA:</label>
                    <input type="number" id="mA">
                </div>
                <div>
                    <label for="mB">mB:</label>
                    <input type="number" id="mB">
                </div>
                <div>
                    <label for="r">r:</label>
                    <input type="number" id="r">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateUniversalGravitation()">Calculate</button>
                </div>
                `;
                break;
            case 'inverse_squares':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Inverse Squares</h2>
                <div>
                    <label for="L">L:</label>
                    <input type="number" id="L">
                </div>
                <div>
                    <label for="d">d:</label>
                    <input type="number" id="d">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateInverseSquares()">Calculate</button>
                </div>
                `;
                break;
            case 'vacuum_frequency':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Vacuum Frequency</h2>
                <div>
                    <label for="lambda_val">λ:</label>
                    <input type="number" id="lambda_val">
                </div>
                <div>
                    <label id="result" class="result">Result:</label>
                </div>
                <button onclick="calculateVacuumFrequency()">Calculate</button>
                </div>
                `;
                break;
            case 'hubble_time':
                contentPlaceholder.innerHTML = `
                <div class="modal-container-user">
                <h2>Calculate Hubble Time CONSTANT</h2>
                <div>
                    <label id="result" class="result">Result: 13968571428.571428</label>
                </div>
                `;
                break;
            
            
    default:
      contentPlaceholder.innerHTML = "<p>No content available.</p>";
  }

  modal.style.display = "block";
}
function closeModal() {
    var modal = document.getElementById("myModal");
    modal.style.display = "none";
  }
  
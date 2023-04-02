# stdlib
import std / [times, strformat, os, parseutils, math, random, strutils]
# nimble
import pkg / [cppstl, arraymancer, shell, unchained]

import garfield, garfield_tools / magboltz_gas

type
  Context = object
    job: int
    date: DateTime # start date
    runNumber: int
    nEvents: int
    dir: string # output path
    degradDir: string # path for Degrad files
    gas: GasMixture
    energy = 10000.0.eV
    length: CentiMeter
    sensor: ptr SensorObj

  DegradHeader = object
    nevt: int
    nclus: int
    nstexc: int
    mcomp: int
    mpair: int
    n4: int
    n5: int
    n6: int
    n7: int
    n8: int
    n9: int
    n10: int
  DegradElement = object
    x: float
    y: float
    z: float
    t: float
    n1: int
    n2: int
    n3: int
  DegradOutput = object
    header: DegradHeader
    data: seq[DegradElement]

# Initialize the pixelmatrix
const Pixelsize = 0.0055
const Pixel = 256 #in one direction
# some other global constants
const DegradOutFile = "DEGRAD.OUT"

randomize(1337)

# forward declare gas ID
proc initContext(job: int,
                 date: DateTime,
                 runNumber: int,
                 nEvents: int,
                 dir, degradDir: string,
                 gas: GasMixture,
                 energy: eV,
                 length: CentiMeter,
                 sensor: ptr SensorObj): Context =
  result = Context(job: job, date: date, runNumber: runNumber,
                   nEvents: nEvents,
                   gas: gas,
                   dir: dir, degradDir: degradDir,
                   energy: energy,
                   length: length,
                   sensor: sensor)


## ##############################
## Generic helpers
## ##############################

proc fixedLength(value: int): string = result = &"{value:06}"

proc formatDate(dt: DateTime): string =
  result = dt.format("yyddMM'_'HH-mm-ss")

proc genOutfile(ctx: Context, event: int): string =
  result = &"run_{fixedLength(ctx.runNumber)}_data_{fixedLength(event)}_{now().formatDate()}"

proc genOutdir(runNumber: int, dt: DateTime): string =
  result = &"Run_{fixedLength(runNumber)}_{dt.formatDate()}"

## ##############################
## Degrad related procs
## ##############################
proc writeDegradFile(fname: string, energy: eV, gas: GasMixture) =
  ## Write an input file for degrad based on photon energy, the gas mixture (two gases and their percentages),
  ## the gas pressure and temperature, the electric field.
  ## It is written into a file with a given filename
  # Get a random seed for degrad
  let seed = rand(1000000)
  echo "[INFO] Degrad seed: ", seed

  # As degrad needs the right amount of spaces fill the file with spaces based on the length of the numbers
  let spaces = @[10, 10,         # row 1
                  5,  5,         # row 2
                 10, 10, 11, 11, # row 3
                 10]             # row 4
  ## We use a template file that has the correct spacing except in places of `$#`
  ## where we will insert the text according to the `spaces` seq
  const tmpl = """
         2         1         3        -1$#$#     2.0000   100.0000
$#$#   99   99   99   99
$#$#     0.000     0.000     0.000     0.000$#$#
$#      0.00      0.00    1    1
    100.00     0.500    0    0    1    1    1    1    1

"""
  let gas1 = gas.gases[0]
  let gas2 = gas.gases[1]
  let data = tmpl % [
    align($seed,                                                     spaces[0]), # row 1
    align(formatFloat(energy.float, ffDecimal, 3),                   spaces[1]), # row 1
    align($(gas1.gID),                                               spaces[2]), # row 2
    align($(gas2.gID),                                               spaces[3]), # row 2
    align(formatFloat(gas1.perc, ffDecimal, 3),                      spaces[4]), # row 3
    align(formatFloat(gas2.perc, ffDecimal, 3),                      spaces[5]), # row 3
    align(formatFloat(gas.temperature.float - 273.15, ffDecimal, 3), spaces[6]), # row 3
    align(formatFloat(gas.pressure.float, ffDecimal, 3),             spaces[7]), # row 3
    align(formatFloat(gas.efield.float, ffDecimal, 2),               spaces[8])] # row 4
  writeFile(fname, data)

#/**
# * Run degrad with a given input file
#*/
proc runDegrad(filename: string) =
  const DegradOutFile2 = "DEGRAD.OUT".toLowerAscii()
  var
    res: string
    err: int
  (res, err) = shellVerbose:
    rm -f ($DegradOutFile)
    rm -f ($DegradOutFile2)
  if err != 0:
    raise newException(IOError, "Could not remove one of the Degrad output files.")
  # manually construct command due to pipe, which is not supported in `shell`
  let cmd = &"./degrad.run < {filename} > degrad.out"
  (res, err) = shellVerbose:
    ($cmd)
  if err != 0:
    raise newException(Exception, "Degrad failed running " & $filename)

proc cleanupDegrad(degrad_dir: string, runNumber, event: int) =
  ## Copy old degrad output file to storage
  let fRun = fixedLength(runNumber)
  let fEv = fixedLength(event)
  let target = &"{degrad_dir}/run_{fRun}_data_{fEv}_{now().formatDate()}.OUT"
  copyFile(DegradOutFile, target)

proc parseDegradOutput(f: string): DegradOutput =
  ## NOTE: If performance was essential here we should instead work with a
  ## mmap'd file, but given that parsing this file is irrelevant for the whole
  ## perf, just read it into a full string for simplicity
  proc parseHeader(s: string): (int, DegradHeader) =
    var idx = 0
    var row = 0
    # we skip all empty space, then parse the integer of the header
    # `fieldPairs` unrolls the loop and we fill each data field of the object
    # one after another
    var header: DegradHeader
    for field, val in fieldPairs(header):
      idx += skipWhile(s, toSkip = {' ', '\t'}, start = idx)
      idx += parseInt(s, val, start = idx)
    result = (idx, header)
  proc parseData(s: string, idx: int): seq[DegradElement] =
    # we skip all empty space, then parse the integer of the header
    # `fieldPairs` unrolls the loop and we fill each data field of the object
    # one after another
    var buf: DegradElement
    var idx = idx # local mutable copy
    while idx < s.len:
      if s[idx] in {'\n', '\l', '\r'}:
        inc idx
      for field, val in fieldPairs(buf):
        idx += skipWhile(s, toSkip = {' ', '\t'}, start = idx)
        when typeof(val) is float:
          idx += parseFloat(s, val, start = idx)
        else:
          idx += parseInt(s, val, start = idx)
      result.add buf
  var sData = readFile(f)
  let (idx, header) = parseHeader(sData)
  let data = parseData(sData, idx+1) # +1 to skip next line!
  result = DegradOutput(header: header,
                        data: data)

proc printDegradOutputInfo(dOut: DegradOutput) =
  echo &"DEGRAD: nevt   = {dOut.header.nevt}"
  echo &"DEGRAD: nclus  = {dOut.header.nclus}"
  echo &"DEGRAD: nstexc = {dOut.header.nstexc}"
  echo &"DEGRAD: mcomp  = {dOut.header.mcomp}"
  echo &"DEGRAD: mpair  = {dOut.header.mpair}"
  echo &"DEGRAD: n4     = {dOut.header.n4}"
  echo &"DEGRAD: n5     = {dOut.header.n5}"
  echo &"DEGRAD: n6     = {dOut.header.n6}"
  echo &"DEGRAD: n7     = {dOut.header.n7}"
  echo &"DEGRAD: n8     = {dOut.header.n8}"
  echo &"DEGRAD: n9     = {dOut.header.n9}"
  echo &"DEGRAD: n10    = {dOut.header.n10}"


## ##############################
## Sampling helpers
## ##############################
import helpers / sampling_helper
func polya*(p: seq[float], x: float): float =
  ## Polya function to fit to TOT histogram / charge in electrons of a
  ## run. This is the actual implementation of the polya distribution.
  ## Parameters:
  ## N     = p[0]    scaling factor
  ## G     = p[1]    gas gain
  ## theta = p[2]    parameter, which describes distribution (?! I guess it makes sens
  ##                 since we take its power and it enters gamma)
  ##
  ## Description goes back to
  ## `Statisitcs of electron avalanches and ultimate resolution of proportional counters`
  ## by Alkhazov, 1970.
  let
    thetaDash = p[2] + 1
    coeff1 = (p[0] / p[1]) * pow((thetaDash), thetaDash) / gamma(thetaDash)
    coeff2 = pow((x / p[1]), p[2]) * exp(-thetaDash * x / p[1])
  result = coeff1 * coeff2

proc initPolyaSampler(): Sampler =
  ## Helper to sample from Pólya
  # For the gas gain the gain is drawn from a polya distribution
  # `TF1 polya = TF1("polya","([0] / [1]) *(((([1]*[1])/([2]*[2]))^(([1]*[1])/([2]*[2]))) /(TMath::Gamma((([1]*[1])/([2]*[2]))))) * ((x /[1])^((([1]*[1])/([2]*[2]))-1)) * exp(-(([1]*[1])/([2]*[2])) *(x / [1]))", 1000, 50000)`
  let fnSample = (
    proc(x: float): float =
      let params = @[4684857630.94, #scaling
                     14076.43, #gain
                     8594.53] #width
      result = polya(params, x)
  )
  # we want to be able to sample between 1k and 50k
  result = sampler(fnSample, 1000.0, 50000.0, num = 5000)

proc makeConfig(dir: string) =
  ## Copies the config files provided in an extra folder into the run folder
  template copyFileIfExist(src, dest: untyped): untyped =
    if existsFile(src):
      copyFile(src, dest)
  #threshold
  copyFileIfExist("config_vorlagen/chip_1_board_0_fec_0_threshold.txt",
                  dir / "chip_1_board_0_fec_0_threshold.txt")
  #matrix
  copyFileIfExist("config_vorlagen/chip_1_board_0_fec_0_matrix.txt",
                  dir / "chip_1_board_0_fec_0_matrix.txt")
  #fsr
  copyFileIfExist("config_vorlagen/chip_1_board_0_fec_0_fsr.txt",
                  dir / "/chip_1_board_0_fec_0_fsr.txt")

proc writeTextToPositionFile(file: string, event: int, energy: eV, angle: Radian,
                                 position: CentiMeter, nclus: int) =
  let data = &"{event}\t{energy.float}\t{angle.float}\t{position.float}\t{nclus}\n"
  var f = open(file, fmAppend)
  f.write(data)
  f.close()

proc getStartPosition(energy: eV, sensor: ptr SensorObj, z_start: CentiMeter): CentiMeter =
  ## Simulate the conversion point of a photon with a given energy and a start point
  ## Use a defined Garfield sensor (with gas mixture, temperature and pressure)
  # Initialize the photon track
  var trackPhoton = TrackHeed.init()
  trackphoton.setSensor(sensor)
  trackphoton.enableElectricField()

  # Initial coordinates of the photon.
  var x0 = 0.0
  var y0 = 0.0
  var z0 = z_start.float
  var t0 = 0.0
  var ne = 0

  # Simulate the photon track. We only care about `ne`
  trackphoton.transportPhoton(x0, y0, z0, t0, energy.float, 0.0, 0.0, -1.0, ne)
  var z_sum = 0.0

  # Iterate over the generated primary electrons to get their starting points
  for i in 0 ..< ne:
    # Variables for the start values of the electrons
    var sx = 0.0
    var sy = 0.0
    var sz = 0.0
    var stime = 0.0
    # Further variables - needed but not meaningful for this application
    var senergy = 0.0
    var sdx = 0.0
    var sdy = 0.0
    var sdz = 0.0
    # Fill the variables for the individual electrons
    trackphoton.getElectron(i, sx, sy, sz, stime, senergy, sdx, sdy, sdz)

    #var drift = AvalancheMC.init()
    #drift.setSensor(sensor)
    #drift.setDistanceSteps(1.μm)
    #drift.driftElectron(sx / 10000.0, sy / 10000.0, sz, stime)


    z_sum += sz
  if ne == 0:
    result = NaN.cm
  else:
    # Get and return the mean starting position of the electrons
    result = (z_sum / ne.float).cm

proc writeEvent(f: string, hits: Tensor[float], number: int) =
  ## Write the output event file
  var data = newStringOfCap(5_000) # only an initial size
  data.add "FEC 0\n"
  data.add "Board 0\n"
  data.add "Chip 1, Hits: {number}\n"
  for x in 0 ..< hits.shape[0]:
    for y in 0 ..< hits.shape[1]:
      if hits[x, y] != 0.0:
        data.add &"{y} {x} {hits[y, x].int}"
  writeFile(f, data)

proc getNextPosition(ctx: var Context, event: var int): CentiMeter =
  ## Retrieve the next position based on the `Sensor`
  while true:
    result = get_start_position(ctx.energy, ctx.sensor, ctx.length)
    if classify(result.float) != fcNaN and result <= ctx.length:
      break
    else:
      inc event

proc propagateElectronGarfield(ctx: var Context, hits: var Tensor[float],
                               number: var int,
                               position: CentiMeter,
                               data: DegradElement,
                               polya: Sampler) =
  ## Performs the propagation of each electron using Garfield++
  # Initialize the electron track for drifting
  #type T = AvalancheMicroscopic
  type T = AvalancheMC
  when T is AvalancheMicroscopic:
    var aval = AvalancheMicroscopic.init() #newAvalancheMicroscopic()
    #AvalancheMC* aval = new AvalancheMC()    # Different simulation tool - not working yet

    aval.setSensor(ctx.sensor)
    #aval->EnablePlotting(&view)  # Debug plot

    # Drift the electron with a given x, y, z start position
    aval.avalancheElectron(data.x / 10000.0, data.y / 10000.0, position.float + (data.z / 10000.0), 0, 0, 0, 0, 0)
    #aval->AvalancheElectron(x / 10000.0, y / 10000.0, position + (z / 10000.), 0)  # Different simulation tool - not working yet
    #aval->DriftElectron(x / 10000.0, y / 10000.0, position + (z / 10000.), 0.0, 0.0, 0.0, 0.0, 0.0) # Different simulation tool - not working yet

    # Get the electron endpoints
    let nElectrons = aval.getNumberOfElectronEndpoints()
    # we only drifted a single electron
    doAssert nElectrons == 1
    var
      x1, y1, z1, t1, e1, x2, y2, z2, t2, e2: float
      status: int
    aval.getElectronEndpoint(0, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, status)
    #aval->GetElectronEndpoint(0, x1, y1, z1, t1, x2, y2, z2, t2, status) # Different simulation tool - not working yet
  else:
    var aval = AvalancheMC.init()    # Different simulation tool - not working yet
    aval.setSensor(ctx.sensor)
    aval.setDistanceSteps(1.μm)
    #aval.enableSignalCalculation()

    #Component* GetComponent(const unsigned int i);

    var
      x0, y0, z0 = 0.005
      ex, ey, ez: float
      medium: ptr Medium
      status0: cint
    ctx.sensor[].electricField(x0, y0, z0, ex, ey, ez, medium, status0)
    echo "\n", medium.isNil, " and status ", status0
    var
      v0, v1, v2: float
    echo medium[].electronVelocity(ex, ey, ez, x0, y0, z0, v0, v1, v2)
    echo ex, " ", ey, " ", ez, " and vel ", v0, " ", v1, " ", v2

    # Drift the electron with a given x, y, z start position
    #aval.driftElectron(data.x / 10000.0, data.y / 10000.0, data.z / 10000.0, 0.0)
    #aval.avalancheElectron(data.x / 10000.0, data.y / 10000.0, position.float + (data.z / 10000.0), 0.0)
    echo data.x, " y ", data.y, " z ", data.z
    aval.driftElectron(data.x / 10000.0, data.y / 10000.0, position.float + (data.z / 10000.0), data.t)
    #aval.avalancheElectron(data.x / 10000.0, data.y / 10000.0, position.float + (data.z / 10000.0), 0, 0, 0, 0, 0)
    #aval->DriftElectron(x / 10000.0, y / 10000.0, position + (z / 10000.), 0.0, 0.0, 0.0, 0.0, 0.0) # Different simulation tool - not working yet

    # Get the electron endpoints
    let nElectrons = aval.getNumberOfElectronEndpoints()
    # we only drifted a single electron
    doAssert nElectrons == 1, "We got " & $nElectrons & " instead of 1!"
    var
      x1, y1, z1, t1, x2, y2, z2, t2: float
      status: int
    aval.getElectronEndpoint(0, x1, y1, z1, t1, x2, y2, z2, t2, status)
    #aval->GetElectronEndpoint(0, x1, y1, z1, t1, x2, y2, z2, t2, status) # Different simulation tool - not working yet


  # Draw a gas gain from the polya distribution
  let amp = polya.sample()
  # Get the pixelcoordinates oh the electron
  let
    posx = floor((x2 + 0.7) / Pixelsize).int
    posy = floor((y2 + 0.7) / Pixelsize).int
  if posx < 0 or posx > 255 or posy < 0 or posy > 255:
    return # nothing to do for this one
  elif hits[posx, posy] == 0: #Counts activated pixels
    inc number
  # Store the hit and the gain in a matrix
  hits[posx, posy] += amp

proc processEvents(ctx: var Context) =
  var event = 0
  var hits = zeros[float]([Pixel, Pixel])
  const angle = 0.0
  for photoelectrons in 0 ..< ctx.nEvents:
    var number = 0 #of activated pixels
    # Clean up the pixelmatrix for a new event
    hits = zeros[float]([Pixel, Pixel])

    echo "[INFO]: photon ", event
    # Generate the filenames for the current event
    let filename = genOutfile(ctx, event)
    let in_file = ctx.degradDir / filename & ".in"

    # Create the degrad in file and run degrad with it
    writeDegradFile(in_file, ctx.energy, ctx.gas)
    runDegrad(in_file)
    echo "[INFO]: Finished degrad"

    # Get a photon conversion point for the event. As an electron might not convert within the detector repeat until a position within the detector is found
    let position = ctx.getNextPosition(event)
    #let position = 0.0.cm
    echo "Final position: ", position

    # Open the degrad output file
    let dOut = parseDegradOutput(DegradOutFile)
    # Read the event parameters from the degrad output file
    printDegradOutputInfo(dOut)

    # Store the truth information of the event in a file
    let photoelectron_file = filename.replace(".txt", "_photoelectrons.txt")
    write_text_to_position_file(photoelectron_file, event, ctx.energy, angle, position, dOut.header.nclus)
    # Iterate over all secondary electrons of the photoelectron track and drift them to the readout
    let polya = initPolyaSampler()
    for iclus in 0 ..< dOut.header.nclus:
      stdout.write("\r[INFO] GARFIELD: Electron: ", iclus + 1, " of ", dOut.header.nclus)
      stdout.flushFile()
      ctx.propagateElectronGarfield(hits, number,
                                    position, dOut.data[iClus],
                                    polya)
    # Close the degrad output file and move it in the runfolder with a new name based on the eventnumber
    cleanupDegrad(ctx.degradDir, ctx.runNumber, event)

    # Write the TOS data file with the zerosupressed x, y and gain data
    let file = ctx.dir / filename & ".txt"
    file.writeEvent(hits, number)
    inc event

proc main(job: int,
          gas1 = "He",
          gas2 = "DME",
          percentage1 = 0.8,
          percentage2 = 0.2,
          temperature = 20.0, # in Celsius!
          pressure = 787.6.Torr,
          eField = 700.0.V•cm⁻¹,
          energy = 10000.0.eV,
          runNumber = 1000,
          nEvents = 10000,
          gasFile = "",
          gasFileDir = "resources/"
         ) =
  # parameters for the simulation - should be arguments in the future
  let gasMixture = initGasMixture((273.15 + temperature).Kelvin,
                                  pressure,
                                  eField,
                                  (gas1, percentage1), (gas2, percentage2))
  let gasFile = readOrGenGasFile(gasMixture, gasFile, gasFileDir)

  # Create a runNumber based on the energy and a job number
  let runNumber = energy.int + job

  # Create the needed folders
  # Get the timestamp - needed for folder- and filenames
  let date = now()
  let dir = genOutdir(runNumber, date)
  let degradDir = dir & "_degrad"
  discard existsOrCreateDir(dir)
  discard existsOrCreateDir(degradDir)
  makeConfig(dir)

  # Make a gas medium.
  let gas = setupMediumMagboltz(gasMixture)

  # Create a cylinder in which the x-rays can convert.
  # Diameter [cm]
  const diameter = 7.8
  # Half-Length of the drfit cylinder [cm] (minus 0.5 cm).
  const length = 2.0.cm
  let tube = SolidTube.init(0.0, 0.0, length.float, 0.5 * diameter, length.float)

  # Combine gas and box to a simple geometry.
  var geo: GeometrySimple
  geo.addSolid(addr tube, gas.get())

  # Make a component with constant electric field.
  var field = ComponentConstant.init()
  field.setGeometry(addr geo)
  field.setElectricField(0.0, 0.0, efield.float)

  # Make a sensor.
  var sensor: SensorObj # = Sensor.new()
  sensor.addComponent(addr field)

  var ctx = initContext(
    job, date, runNumber, nEvents,
    dir, degradDir,
    gasMixture,
    energy,
    length,
    addr sensor)
  processEvents(ctx)

when isMainModule:
  import unchained / cligenParseUnits
  import cligen
  dispatch(main, help = {
    "job" : "The job ID for this process",
    "nEvents" : "Number of events to simulate",
    "gas1" : "The main gas to use",
    "gas2" : "The secondary gas to use",
    "percentage1" : "Percentage of gas 1 (as a fraction of 1)",
    "percentage2" : "Percentage of gas 2 (as a fraction of 1)",
    "temperature" : "Temperature in Celsius of the gas",
    "pressure" : "Pressure of the gas in Torr",
    "eField" : "Electric field applied to the chamber",
    "energy" : "Energy of the X-rays to simulate in eV",
    "runNumber" : "Run number to use for this job as a start (overwritten)",
    "gasFile" : """The gas file to use. Should match the gas 1 and gas 2 (type & percentages).
  If none given will generate a gas file.""",
    "gasFileDir" : """The location where gas files are stored. If no gas file given will first
  attempt to read a corresponding gas file from there"""
  })

import os, sys, re, json, glob, shutil, fcntl, datetime
from contextlib import contextmanager

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from parsing import parse_metadata
from parameters import MetadataCache

RUN_PREFIX = "run"
RUN_NAME_REGEX = re.compile(rf"^{RUN_PREFIX}(\d+)$", re.IGNORECASE)

class Run:
    '''
    An analysis run holds the segregation statistics produced by psQTL_proc.py (and any
    outputs from psQTL_post.py) for one particular grouping of samples. Runs are stored
    as numbered subdirectories (run1, run2, ...) of the working directory's 'runs' folder,
    with each containing a copy of the metadata file that defines its groups.
    '''
    infoFileName = "run.json"
    metadataFileNames = ["metadata.tsv", "metadata.tsv.gz"]

    def __init__(self, runDir):
        self.runDir = os.path.abspath(runDir)
        self.number = run_number_from_name(os.path.basename(self.runDir))
        if self.number == None:
            raise ValueError(f"'{self.runDir}' is not named like an analysis run directory (e.g., 'run1')")
        self._metadataDict = None
        self._info = None

    @property
    def name(self):
        return f"{RUN_PREFIX}{self.number}"

    @property
    def metadataFile(self):
        for fileName in Run.metadataFileNames:
            metadataFile = os.path.join(self.runDir, fileName)
            if os.path.isfile(metadataFile):
                return metadataFile
        return None

    @property
    def metadataDict(self):
        if self._metadataDict == None:
            if self.metadataFile == None:
                raise FileNotFoundError(f"Analysis run '{self.runDir}' is missing its metadata file!")
            self._metadataDict = parse_metadata(self.metadataFile)
        return self._metadataDict

    @property
    def infoFile(self):
        return os.path.join(self.runDir, Run.infoFileName)

    @property
    def info(self):
        if self._info == None:
            self._info = {}
            if os.path.isfile(self.infoFile):
                with open(self.infoFile, "r") as fileIn:
                    try:
                        self._info = json.load(fileIn)
                    except json.JSONDecodeError:
                        pass # run details are informational only; the metadata file is what matters
        return self._info

    @property
    def metadataSource(self):
        return self.info.get("metadataSource", None)

    @property
    def created(self):
        return self.info.get("created", None)

    @property
    def migrated(self):
        return self.info.get("migrated", False)

    @property
    def sourceChanged(self):
        '''
        Returns True if the metadata file this run was created from now indicates different
        groups, False if it is unchanged, or None if the source file no longer exists.
        '''
        if self.metadataSource == None or not os.path.isfile(self.metadataSource):
            return None
        try:
            return not self.matches(parse_metadata(self.metadataSource))
        except ValueError:
            return True

    def matches(self, metadataDict):
        "parse_metadata() returns sorted sample lists, so equality is order-independent"
        return self.metadataDict["group1"] == metadataDict["group1"] and \
               self.metadataDict["group2"] == metadataDict["group2"]

    def describe(self):
        '''
        Returns a one-line description of this run suitable for listing runs to the user.
        '''
        description = (f"{self.name}: group1 n={len(self.metadataDict['group1'])}, " +
                       f"group2 n={len(self.metadataDict['group2'])}")
        if self.metadataSource != None:
            description += f"; metadata from '{self.metadataSource}'"
            if self.sourceChanged == True:
                description += " (since modified)"
        if self.migrated:
            description += "; migrated from an earlier psQTL version"
        elif self.created != None:
            description += f"; created {self.created}"
        return description

def run_number_from_name(name):
    match = RUN_NAME_REGEX.match(name)
    return int(match.group(1)) if match != None else None

def parse_run_id(value):
    '''
    Parameters:
        value -- a string or integer provided by the user to identify a run e.g., '2' or 'run2'
    Returns:
        runNumber -- an integer indicating the run number
    '''
    value = str(value).strip()
    runNumber = int(value) if value.isdigit() else run_number_from_name(value)
    if runNumber == None or runNumber < 1:
        raise ValueError(f"--run value '{value}' is not recognised; specify a run number like '2' or 'run2'")
    return runNumber

def list_runs(workingDirectory):
    '''
    Returns:
        runs -- a list of Run objects sorted by run number; directories without a metadata
                file are ignored since they cannot be a complete run
    '''
    runsDir = os.path.join(workingDirectory, "runs")
    if not os.path.isdir(runsDir):
        return []
    runs = []
    for name in os.listdir(runsDir):
        runDir = os.path.join(runsDir, name)
        if os.path.isdir(runDir) and run_number_from_name(name) != None:
            run = Run(runDir)
            if run.metadataFile != None:
                runs.append(run)
    return sorted(runs, key=lambda run: run.number)

def describe_runs(runs):
    if runs == []:
        return "No analysis runs exist yet."
    return "Existing analysis runs:\n" + "\n".join([f"    {run.describe()}" for run in runs])

def get_run(workingDirectory, runNumber):
    runs = list_runs(workingDirectory)
    for run in runs:
        if run.number == runNumber:
            return run
    raise FileNotFoundError(f"Analysis run '{RUN_PREFIX}{runNumber}' does not exist in '{workingDirectory}'.\n" +
                            describe_runs(runs))

def find_matching_run(workingDirectory, metadataDict):
    for run in list_runs(workingDirectory):
        if run.matches(metadataDict):
            return run
    return None

@contextmanager
def runs_lock(workingDirectory):
    '''
    Serialises run creation so that, for example, simultaneously launched 'psQTL_proc.py ed'
    and 'psQTL_proc.py splsda' commands with the same new metadata share one run.
    '''
    runsDir = os.path.join(workingDirectory, "runs")
    os.makedirs(runsDir, exist_ok=True)
    with open(os.path.join(runsDir, ".lock"), "a") as lockFile:
        fcntl.flock(lockFile, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lockFile, fcntl.LOCK_UN)

def _create_run(workingDirectory, metadataDict, sourceFile=None, migrated=False):
    '''
    Creates the next numbered run directory. Must be called within runs_lock().

    Parameters:
        workingDirectory -- a string indicating the psQTL working directory
        metadataDict -- a dictionary of groups as returned by parse_metadata()
        sourceFile -- (OPTIONAL) a string indicating the metadata file the groups were
                      parsed from; it is copied verbatim so that sPLS-DA reads exactly
                      what the user provided. If None, a metadata file is written from
                      metadataDict instead.
        migrated -- (OPTIONAL) a boolean indicating whether this run is being created to
                    hold results from an earlier version of psQTL
    '''
    runsDir = os.path.join(workingDirectory, "runs")
    existingNumbers = [ run_number_from_name(name) for name in os.listdir(runsDir) ]
    runNumber = max([ num for num in existingNumbers if num != None ], default=0) + 1
    runDir = os.path.join(runsDir, f"{RUN_PREFIX}{runNumber}")
    os.mkdir(runDir)

    # Store run details
    info = {
        "metadataSource": os.path.abspath(sourceFile) if sourceFile != None else None,
        "created": datetime.datetime.now().isoformat(sep=" ", timespec="seconds"),
        "migrated": migrated
    }
    with open(os.path.join(runDir, Run.infoFileName), "w") as fileOut:
        json.dump(info, fileOut, indent=4)

    # Store metadata last, since its presence marks the run as complete
    metadataFileName = Run.metadataFileNames[1] if sourceFile != None and sourceFile.endswith(".gz") \
                       else Run.metadataFileNames[0]
    tmpFile = os.path.join(runDir, f".{metadataFileName}.tmp")
    if sourceFile != None:
        shutil.copyfile(sourceFile, tmpFile)
    else:
        with open(tmpFile, "w") as fileOut:
            for group in ["group1", "group2"]:
                for sample in metadataDict[group]:
                    fileOut.write(f"{sample}\t{group}\n")
    os.replace(tmpFile, os.path.join(runDir, metadataFileName))
    return Run(runDir)

def find_or_create_run(workingDirectory, metadataFile):
    '''
    Returns:
        run -- a Run object whose groups match those in metadataFile
        created -- a boolean indicating whether the run was newly created
    '''
    metadataDict = parse_metadata(metadataFile)
    with runs_lock(workingDirectory):
        run = find_matching_run(workingDirectory, metadataDict)
        if run != None:
            return run, False
        return _create_run(workingDirectory, metadataDict, sourceFile=metadataFile), True

def find_legacy_results(workingDirectory):
    '''
    psQTL versions prior to analysis runs wrote ED files into the working directory
    itself and sPLS-DA files into its 'splsda' subdirectory.

    Returns:
        legacyFiles -- a sorted list of ED (and associated .ok and .pkl) files
        legacySplsdaDir -- a string indicating the legacy 'splsda' directory, or None
    '''
    legacyFiles = set()
    for pattern in ["psQTL_call.*_ed.*", "psQTL_depth.ed.*"]:
        legacyFiles.update([ f for f in glob.glob(os.path.join(workingDirectory, pattern)) if os.path.isfile(f) ])
    legacySplsdaDir = os.path.join(workingDirectory, "splsda")
    return sorted(legacyFiles), legacySplsdaDir if os.path.isdir(legacySplsdaDir) else None

def migrate_legacy_results(workingDirectory, defaultMetadataFile=None):
    '''
    Moves results from an earlier psQTL version into an analysis run. The metadata used to
    generate them is assumed to be the metadata given to psQTL_prep.py initialise, falling
    back to the groups recorded in the metadata cache if that file no longer exists.

    Returns:
        run -- the Run object that results were moved into, or None if there was nothing to migrate
    '''
    legacyFiles, legacySplsdaDir = find_legacy_results(workingDirectory)
    if legacyFiles == [] and legacySplsdaDir == None:
        return None

    # Determine the groups that legacy results were generated with
    sourceFile, metadataDict = None, None
    if defaultMetadataFile != None and os.path.isfile(defaultMetadataFile):
        sourceFile = defaultMetadataFile
        metadataDict = parse_metadata(sourceFile)
    else:
        try:
            metadataCache = MetadataCache(workingDirectory)
            metadataCache.load()
            if metadataCache.group1 != None and metadataCache.group2 != None:
                metadataDict = { "group1": sorted(metadataCache.group1), "group2": sorted(metadataCache.group2) }
        except FileNotFoundError:
            pass
    if metadataDict == None:
        raise ValueError(f"Results from an earlier version of psQTL were found in '{workingDirectory}', but the " +
                         "metadata they were generated with could not be determined. Move or delete the " +
                         "'psQTL_call*' and 'psQTL_depth*' files and the 'splsda' directory, then re-run psQTL_proc.py.")

    with runs_lock(workingDirectory):
        run = find_matching_run(workingDirectory, metadataDict)
        if run == None:
            run = _create_run(workingDirectory, metadataDict, sourceFile=sourceFile, migrated=True)

        # Check for conflicts before moving anything
        moves = [ [f, os.path.join(run.runDir, os.path.basename(f))] for f in legacyFiles ]
        if legacySplsdaDir != None:
            moves.append([legacySplsdaDir, os.path.join(run.runDir, "splsda")])
        conflicts = [ dst for _, dst in moves if os.path.exists(dst) ]
        if conflicts != []:
            raise FileExistsError(f"Results from an earlier version of psQTL were found in '{workingDirectory}' but " +
                                  f"could not be moved into '{run.runDir}' since it already contains: " +
                                  f"{', '.join(conflicts)}. Resolve this manually before re-running.")
        for src, dst in moves:
            shutil.move(src, dst)

    print(f"# Moved results from an earlier version of psQTL into '{run.runDir}'; these were " +
          f"generated with metadata from '{sourceFile if sourceFile != None else 'the metadata cache'}'")
    return run

def resolve_run(locations, runID=None, metadataFile=None, defaultMetadataFile=None, create=False):
    '''
    Locates the analysis run that psQTL_proc.py or psQTL_post.py should operate within,
    and sets it as the run directory of the provided Locations object.

    Precedence is given to runID, then metadataFile, then defaultMetadataFile (i.e., the
    metadata given to psQTL_prep.py initialise), and finally the sole existing run.

    Parameters:
        locations -- a Locations object for the working directory
        runID -- (OPTIONAL) a string from --run identifying an existing run e.g., '2' or 'run2'
        metadataFile -- (OPTIONAL) a string from --meta indicating a metadata file
        defaultMetadataFile -- (OPTIONAL) a string indicating the cached metadata file
        create -- (OPTIONAL) a boolean indicating whether a new run should be created if no
                  existing run matches the metadata
    Returns:
        run -- the resolved Run object
    '''
    workingDirectory = locations.workingDirectory
    if runID != None and metadataFile != None:
        raise ValueError("--run and --meta cannot be specified together; provide one or the other")

    migrate_legacy_results(workingDirectory, defaultMetadataFile)

    # Identify which metadata, if any, should be used to locate the run
    if runID == None and metadataFile == None and defaultMetadataFile != None:
        if os.path.isfile(defaultMetadataFile):
            metadataFile = defaultMetadataFile
        else:
            print(f"# Note: metadata file '{defaultMetadataFile}' given to psQTL_prep.py initialise no " +
                  "longer exists; looking for an existing analysis run instead")
    if metadataFile != None and not os.path.isfile(metadataFile):
        raise FileNotFoundError(f"--meta file '{metadataFile}' does not exist or is not a file!")

    # Locate the run
    created = False
    if runID != None:
        run = get_run(workingDirectory, parse_run_id(runID))
    elif metadataFile != None:
        if create:
            run, created = find_or_create_run(workingDirectory, metadataFile)
        else:
            run = find_matching_run(workingDirectory, parse_metadata(metadataFile))
            if run == None:
                raise FileNotFoundError(f"No analysis run matches the groups in metadata file '{metadataFile}'; " +
                                        "run psQTL_proc.py with this metadata first.\n" +
                                        describe_runs(list_runs(workingDirectory)))
    else:
        runs = list_runs(workingDirectory)
        if len(runs) == 1:
            run = runs[0]
        elif len(runs) == 0:
            raise ValueError("No metadata has been provided; specify a metadata file with --meta (or with " +
                             "'psQTL_prep.py initialise --meta') to indicate how samples are grouped.")
        else:
            raise ValueError("Multiple analysis runs exist; specify which one to use with --run or --meta.\n" +
                             describe_runs(runs))

    print(f"# {'Created' if created else 'Using'} analysis run {run.describe()}")
    locations.runDir = run.runDir
    return run

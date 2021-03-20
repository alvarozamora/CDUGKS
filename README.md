<h1>Massively Parallel Coupled Discrete Unified Gas Kinetic Scheme</h1>

Welcome to the MP-CDUGKS github repository. MP-CDUGKS is written in the [Regent](https://regent-lang.org) language, and uses the [Legion Runtime System](https://github.com/StanfordLegion/legion). I recommend using the `control_replication` branch, which as of writing this has better one-node performance for this code with the `-dm:exact` runtime flag.

Refer to the Legion repository for instructions on how to build the runtime system.

<h2>Quick Start Guide </h2>
There are O(~10) test problems currently implemented, with `testProblem` id

1) Sod Shock Tube
2) Kelvin-Helmholtz (Sharp Interface, wavelength = 1)
3) Uniform Shearing Box
4) Ramped Kelvin-Helmholtz
5) Nonuniform Shearing Box
6) Maxwellian Relaxation (Nonstandard problem: Two blobs moving toward each other. Sine Wave Collapse is a better continuous/periodic version of this problem)
7) Cloud Crushing (Blob) Test. Note: No elongated box is used.
8) Thermoacoustic Wave
9) Gresho Vortex

To run one of these problems, run `path/to/regent/executable/regent.py Main.rg -p testProblem -c <subregions> -ll:cpu <cores/node> -ll:csize <mem/node>`. If using the `control_replication` branch, also add the `-dm:exact` flag, which instructs the default mapper to map exact regions to cores when only using one node.

The conserved variables will be output at every timestep to the relative `Data/` path.

<h2>Adding Test Problems</h2>
To add test problems, you will first need to set simulation parameters in the task `TestProblem`. Then, you will need to specify the initial conditions in `InitializeW`. For Non-Maxwellian initializations, you will need to modify `InitializeGrid`.

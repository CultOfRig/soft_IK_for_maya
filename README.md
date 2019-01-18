# soft_IK_for_maya
A simple plugin providing the source for node(s) solving V-Chain IK hierarchies with soft-solution options

This is being worked on for the January 2019 Patreon reward of Cult of Rig.

Source is available to all under BSD 3 Clause, and the binaries for Maya 2018 are provided in the dist directory.
Most people should be able to use the AVX2 binary, but if you have an old CPU (pre-Haswell for intel or pre-Ryzen for AMD) you might have to use the SSE 4.2 (SSE2) version to avoid crashes.

A demo scene with two set ups showing hierarchical and flat modes are available in the demo directory.

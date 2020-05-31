# soft_IK_for_maya
A simple plugin providing the source for node(s) solving V-Chain IK hierarchies with soft-solution options
This is an accessory to January 2019's Patreon reward of Cult of Rig, but the nodes or code themselves could be generally useful.

This is offered for its pedagogic value.
It's not a bad option to use in production, since it's not exactly raining closed-form solvers even for simple cases like a vChain, but it's written and laid out to teach about the solver, not necessarily to be the fastest or most optimal.

If you want to try and speed it up, please, feel free, but do know that the solution part of the node cost is already less than 20% of the evaluation cost, most of the cost in an actual rig will be handle tax and Maya's internal orchestration.

Source is available under BSD 3 Clause, and the binaries for Maya 2018 are provided in the dist directory.
Most people should be able to use the AVX2 binary, but if you have an old CPU (pre-Haswell for intel or pre-Ryzen for AMD) you might have to use the SSE 4.2 (SSE2) version to avoid crashes.

A demo scene is available in the demo directory with two set-up's showing hierarchical and flat solution implementation.

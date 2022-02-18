// Basic vs Advanced algorithms, not utilities

// Basic algorithm:  G can be modified (input/output).  Properties computed
// if needed.  Sometimes fewer inputs than Advanced (parameters computed
// automically, or with defaults).

LAGraph_WhatEver    // is a Basic method

LAGraph_BreadthFirstSearch          // so this is Basic

// Advanced algorithm:  G is input-only.  Error if properties missing.
// sometimes has more input parameters.

LAGraph_BreadthFirstSearch_         // append an underscore;
                                    // I think this is a bit odd

LAGraph_BreadthFirstSearch_Z        // append an underscore and _something;
                                    // "Z" would be consistent across all algo:
LAGraph_BreadthFirstSearch_Adv      // as in "Advanced".  Too wordy

LAGraph__BreadthFirstSearch         // two underscores in the front,
LAGraph_BreadthFirstSearch__xxx     // or at the end.
                                    // LAGraph__WhatEver_MoreStuff is advanced
                                    // Hard to notice.

LAGraph_BreadthFirstSearch_xxx      // xxx is some algo-dependent suffix.
                                    // Too dependent on the particular algo.

LAGraph_BreadthFirstSearch_pushpull // "pushpull" is lower case
                                    // hard to notic

LAGr_BreadthFirstSearch             // different prefix.  I like this best.

LG_BreadthFirstSearch               // LG_* is currently an internal prefix for
                                    // non-user-callable methods.  We could
                                    // revise that, and make LG_* used for
                                    // Advanced methods.


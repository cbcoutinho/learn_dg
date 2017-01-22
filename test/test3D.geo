L = 1;
dx1 = 0.2*L;

// Create points on x-z plane
Point(1) = {-0.5*L, -0.5*L, -0.5*L, dx1};
Point(2) = { 0.5*L, -0.5*L, -0.5*L, dx1};
Point(3) = { 0.5*L, -0.5*L,  0.5*L, dx1};
Point(4) = {-0.5*L, -0.5*L,  0.5*L, dx1};

// Create lines connecting points
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Create loop composing of lines 1-4
ll = newll;
Line Loop(ll) = {1, 2, 3, 4};

// Create surface `sc` from loop `ll`
sc = news;
Plane Surface(sc) = {ll};

// Extrude `sc` along y axis to create box
ex = Extrude {0, 1*L, 0} { Surface{sc}; };

// Use gmsh numbering
//Physical Surface(sc)      = {sc};     // front
//Physical Surface(ex[0])    = {ex[0]};  // back
//Physical Surface(ex[2])  = {ex[2]};  // bottom
//Physical Surface(ex[3])   = {ex[3]};  // right
//Physical Surface(ex[4])     = {ex[4]};  // top
//Physical Surface(ex[5])    = {ex[5]};  // left
//Physical Volume(ex[1]) = {ex[1]};

// use 0-N numbering with named groups
//Physical Surface('front')      = {sc};     // front
//Physical Surface('back')    = {ex[0]};  // back
//Physical Surface('bottom')  = {ex[2]};  // bottom
//Physical Surface('right')   = {ex[3]};  // right
//Physical Surface('top')     = {ex[4]};  // top
//Physical Surface('left')    = {ex[5]};  // left
//Physical Volume('volume') = {ex[1]};

// use both gmsh and 0-N numbering
//Physical Surface('front',sc)      = {sc};     // front
//Physical Surface('back',ex[0])    = {ex[0]};  // back
//Physical Surface('bottom',ex[2])  = {ex[2]};  // bottom
//Physical Surface('right',ex[3])   = {ex[3]};  // right
//Physical Surface('top',ex[4])     = {ex[4]};  // top
//Physical Surface('left',ex[5])    = {ex[5]};  // left
//Physical Volume('volume',ex[1]) = {ex[1]};

// Use 0-N numbering directly
Physical Surface(1)      = {sc};     // front
Physical Surface(2)    = {ex[0]};  // back
Physical Surface(3)  = {ex[2]};  // bottom
Physical Surface(4)   = {ex[3]};  // right
Physical Surface(5)     = {ex[4]};  // top
Physical Surface(6)    = {ex[5]};  // left
Physical Volume(7) = {ex[1]};

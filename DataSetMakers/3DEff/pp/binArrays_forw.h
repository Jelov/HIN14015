static const double yarray[]        = {1.6, 2.4};
static const double ptarray[]       = {3.0, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 11.0, 13.0, 16.0, 30.0};
static const double ctauarray[]     = {0, 0.3, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};
static const double ctauforwarray[] = {0, 0.4, 0.8, 1.2, 1.6, 3.0, 4.0, 5.0, 7.0, 10.0};
static const int centarray[]        = {0, 40};

static const int nbinsy = sizeof(yarray)/sizeof(double);
static const int nbinspt = sizeof(ptarray)/sizeof(double);
static const int nbinsctau = sizeof(ctauarray)/sizeof(double);
static const int nbinsctauforw = sizeof(ctauforwarray)/sizeof(double);
static const int nbinscent = sizeof(centarray)/sizeof(int);

static const double resolmin = -500;
static const double resolmax = 500;
static const int nbinsresol = 125;


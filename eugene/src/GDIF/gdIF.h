extern "C" void InitPNG(int x,int y, int From, int To,int olap, int ILen,char *name);
extern "C" void PlotBarF(unsigned int nuc, signed char phase, REAL pos, REAL width, int col);
extern "C" void PlotBarI(unsigned int nuc, signed char phase, REAL pos, int width, int col);
extern "C" void PlotLine(unsigned int nuc1, unsigned int nuc2, 
			 signed char phase1, signed char phase2, REAL pos1, REAL pos2, int col);
extern "C" void ClosePNG();
extern "C" void OutputHTMLFileNames();

#include "errlog.h"
#include <fstream>
#include <time.h>

Errlog::Errlog(void)
{
}

Errlog::~Errlog(void)
{
}

void Errlog::errlog(const char* errinfor)
{
	time_t nowtime;
    std::ofstream fout;
    fout.open("./Errinfor.log", std::ofstream::out | std::ofstream::app);
	if(fout.is_open())
	{
		nowtime = time(NULL);
        fout << "->" << ctime(&nowtime) << errinfor << std::endl;
	}
	fout.close();
}

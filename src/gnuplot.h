/**
 *  Copyright 2016 by Alexandre Aguiar
 *
 * This file is part of Open Electromagnetic Transient Program - OEMTP.
 * 
 * OEMTP is free software: you can redistribute 
 * it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation, either 
 * version 3 of the License, or (at your option) any later version.
 * 
 * Some open source application is distributed in the hope that it will 
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @license GPL-3.0+ <http://spdx.org/licenses/GPL-3.0+>
 */

#ifndef GNUPLOT_H
#define GNUPLOT_H
#ifdef  GNUPLOT
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <list>
#define _PCLOSE pclose
#define _POPEN  popen

using std::string;
using std::vector;

template<class T>
class gnuplot
{
    public:
        gnuplot(string outputName, string graphT, string outputFolder);
        virtual ~gnuplot();
        bool Join(string xlabel, string ylabel, string legend, vector<T>& ydata, vector<T>& xdata, int breakPlot);
        bool run(int pointInterval, bool colour);

    private:
        string folder;
        string outputFile;
        string graphType;
        vector<string> xlabels;
        vector<string> ylabels;
        vector<string> legends;
        vector<int> breakPlots;
        vector< vector<T>* > dataVectors;
        vector< vector<T> > timeVectors;
        FILE *p;
};

template <class T> gnuplot<T>::gnuplot(string outputName, string graphT, string outputFolder){
    outputFile = outputName;
    graphType  = graphT;
    folder     = outputFolder;
    mkdir(folder.c_str(),S_IRUSR|S_IWUSR|S_IXUSR);
}

template <class T> gnuplot<T>::~gnuplot(){}

template <class T> bool gnuplot<T>::Join(string xlabel,string ylabel, string legend, vector<T>& ydata, vector<T>& xdata, int breakPlot){
	xlabels.push_back(xlabel);
	ylabels.push_back(ylabel);
	legends.push_back(legend);
	dataVectors.push_back(&ydata);
	timeVectors.push_back(xdata);
	breakPlots.push_back(breakPlot);
 	//cout<<breakPlot<<endl;
    return 0;
}

template <class T> bool gnuplot<T>::run(int pointInterval, bool colour){

	int nGraphs = 0, linesPerPlot, j, z, i, key;
	std::list<int> uniquePlots(breakPlots.begin(), breakPlots.end()), indexes, aux_index;
	uniquePlots.unique();
	nGraphs = uniquePlots.size();
	vector<string> command, extension;
	command.push_back("wxt size 1200,800");
	command.push_back("png size 1200,800");
	command.push_back("pdfcairo");
	extension.push_back("png");
	extension.push_back("png");
	extension.push_back("pdf");
	for(unsigned k = 0; k < command.size(); k++){
		uniquePlots.assign(breakPlots.begin(), breakPlots.end());
		uniquePlots.unique();
		for(int graf = 0;graf<nGraphs;graf++){
			linesPerPlot = 0;
			key = uniquePlots.back(); uniquePlots.pop_back();
			for(int j = 0; j< breakPlots.size(); j++){
				if (breakPlots[j] == key){
					linesPerPlot++;
					indexes.push_front(j);					
				}
			}
			aux_index = indexes;
			p = _POPEN("gnuplot -persist","w");
			fprintf(p, "set term %s \n",command[k].c_str());
			fprintf(p, "set output \"%s/%s_%s.%s\" \n",folder.c_str(),outputFile.c_str(),legends[indexes.back()].c_str(),extension[k].c_str());
			if(graphType == "logxy")
				fprintf(p, "set logscale xy \n");
			if(graphType == "logx")
				fprintf(p, "set logscale x \n");
			fprintf(p, "set multiplot layout 1,1 \n"); 
			fprintf(p, "set grid \n"); 
			i = aux_index.back(); aux_index.pop_back();
			fprintf(p, "set xlabel \"%s\" \n", xlabels[i].c_str());
			fprintf(p, "set ylabel \"%s\" \n", ylabels[i].c_str());
			fprintf(p, "plot '-' using 1:2 title \"%s\" with lp pi %d",legends[i].c_str(),pointInterval);
			if(linesPerPlot > 1)
				for(int lPP = 1; lPP<(linesPerPlot);lPP++){
					i = aux_index.back(); aux_index.pop_back();
					fprintf(p, ", '-' using 1:2 title \"%s\" with lp pi %d",legends[i].c_str(),pointInterval);
				}
			fprintf(p,"\n");
			for(z = 0; z < linesPerPlot; z++){
				i = indexes.back(); indexes.pop_back();
				for(unsigned j = 0; j<timeVectors[i].size(); j++)
					fprintf(p, "%lf %lf\n", timeVectors[i][j], (*dataVectors[i])[j]);
				fprintf(p, "e\n");
			}
			fprintf(p,"exit \n");
			_PCLOSE(p);
		}
	}

	uniquePlots.assign(breakPlots.begin(), breakPlots.end());
	uniquePlots.unique();
    	for(int graf = 0;graf<nGraphs;graf++){
		linesPerPlot = 0;
		key = uniquePlots.back(); uniquePlots.pop_back();
		for(int j = 0; j< breakPlots.size(); j++){
			if (breakPlots[j] == key){
				linesPerPlot++;
				indexes.push_front(j);					
			}
		}
		aux_index = indexes;
		p = _POPEN("gnuplot -persist","w");
		if(graphType=="logxy")
			fprintf(p, "set logscale xy \n");
		if(graphType=="logx")
            	fprintf(p, "set logscale x \n");
		fprintf(p,"set term postscript eps \n");
		if (colour)
			fprintf(p, "set term postscript eps enhanced color\n");
		fprintf(p,"set output \"%s/%s_%s.eps\" \n",folder.c_str(),outputFile.c_str(),legends[indexes.back()].c_str());
		fprintf(p, "set multiplot layout 1,1 \n");  
		fprintf(p, "set grid \n");  
		i = aux_index.back(); aux_index.pop_back();
		fprintf(p, "set xlabel \"%s\" \n", xlabels[i].c_str());
		fprintf(p, "set ylabel \"%s\" \n", ylabels[i].c_str());
		fprintf(p, "plot '-' using 1:2 title \"%s\" with lp pi %d",legends[i].c_str(),pointInterval);
		if(linesPerPlot > 1)
			for(int lPP = 1; lPP<(linesPerPlot);lPP++){
				i = aux_index.back(); aux_index.pop_back();
				fprintf(p, ", '-' using 1:2 title \"%s\" with lp pi %d",legends[i].c_str(),pointInterval);
			}
		fprintf(p,"\n");
		for(z = 0; z < linesPerPlot; z++){
			i = indexes.back(); indexes.pop_back();
			for(unsigned j = 0; j<timeVectors[i].size(); j++)
				fprintf(p, "%lf %lf\n", timeVectors[i][j], (*dataVectors[i])[j]);
				fprintf(p, "e\n");
		}

			fprintf(p,"exit \n");
			_PCLOSE(p);
	}
 
return 0;
}

#endif // GNUPLOT_H
#endif // GNUPLOT

/*template <class T> bool gnuplot<T>::run(int pointInterval, bool colour){

	vector<string> command, extension;
	command.push_back("wxt size 1200,800");
	command.push_back("png size 1200,800");
	command.push_back("pdfcairo");
	extension.push_back("png");
	extension.push_back("png");
	extension.push_back("pdf");
	for(unsigned k = 0; k < command.size(); k++)
		for(unsigned i = 0; i < dataVectors.size(); i++){
			p = _POPEN("gnuplot -persist","w");
			if(graphType=="logxy")
				fprintf(p, "set logscale xy \n");
			if(graphType=="logx")
				fprintf(p, "set logscale x \n");
			fprintf(p, "set term %s \n",command[k].c_str());
			fprintf(p, "set output \"%s/%s_%d.%s\" \n",folder.c_str(),outputFile.c_str(),i,extension[k].c_str());
			fprintf(p, "set xlabel \"%s\" \n",xlabels[i].c_str());
			fprintf(p, "set ylabel \"%s\" \n",ylabels[i].c_str());
			if (colour){
				fprintf(p, "plot '-' using 1:2 title \"%s\" with lp pi %d\n",legends[i].c_str(),pointInterval);
			}else{
				fprintf(p, "plot '-' using 1:2 title \"%s\" with lp pi %d lc \"black\" \n",legends[i].c_str(),pointInterval);
			}
			for(unsigned j=0;j<(*dataVectors[i]).size();j++){
				fprintf(p,"%lf %lf\n",timeVectors[i][j],(*dataVectors[i])[j]);
			}
			fprintf(p,"e\n");
			fprintf(p,"exit \n");
			_PCLOSE(p);
		}

	for(unsigned i=0;i<dataVectors.size();i++){
		p = _POPEN("gnuplot -persist","w");
            if(graphType== "logxy")
                fprintf(p, "set logscale xy \n");
            if(graphType== "logx")
                fprintf(p, "set logscale x \n");
            fprintf(p, "set term postscript eps \n");
            if (colour)
                fprintf(p, "set term postscript eps enhanced color\n");
            fprintf(p, "set output \"%s/%s_%d.eps\" \n",folder.c_str(),outputFile.c_str(),i);
            fprintf(p, "set xlabel \"%s\" \n",xlabels[i].c_str());
            fprintf(p, "set ylabel \"%s\" \n",ylabels[i].c_str());
            fprintf(p, "plot '-' using 1:2 title \"%s\" with lp pi %d\n",legends[i].c_str(),pointInterval);
            for(unsigned j=0;j<(*dataVectors[i]).size();j++){
                fprintf(p, "%lf %lf\n", timeVectors[i][j], (*dataVectors[i])[j]);
            }
            fprintf(p, "e\n");
            fprintf(p, "exit \n");
            _PCLOSE(p);
	}
	return 0;
}
*/

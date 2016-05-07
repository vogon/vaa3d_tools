/* consensus_skeleton_func.cpp
 * a plugin to merge multiple neurons by generating a consensus skeleton
 * 2012-05-02 : by Yinan Wan
 */

#include <v3d_interface.h>
#include "basic_surf_objs.h"
#include "v3d_message.h"
#include "consensus_skeleton_func.h"
#include "consensus_skeleton.h"
#include "median_swc.h"
#include "dark_pruning.h"
#include <vector>
#include <iostream>
using namespace std;

const QString title = QObject::tr("Consensus Skeleton");

int consensus_swc_menu(V3DPluginCallback2 &callback, QWidget *parent)
{
	QString fileOpenName;
	fileOpenName = QFileDialog::getOpenFileName(0, QObject::tr("Open File"),
			"",
			QObject::tr("Supported file (*.ano)"
				";;Neuron structure	(*.ano)"
				));
	if(fileOpenName.isEmpty()) 
		return -1;
		
	P_ObjectFileType linker_object;
	if (!loadAnoFile(fileOpenName,linker_object))
	{
		fprintf(stderr,"Error in reading the linker file.\n");
		return -1;
	}
	
	QStringList nameList = linker_object.swc_file_list;
	V3DLONG neuronNum = nameList.size();
	vector<NeuronTree> nt_list;
	
	for (V3DLONG i=0;i<neuronNum;i++)
	{
		NeuronTree tmp = readSWC_file(nameList.at(i));
		nt_list.push_back(tmp);
	}

	bool ok;
    //int method_code = QInputDialog::getInt(parent, "method", "FIX ME: 1:direct vote; 2: maximum spanning tree", 0, 1, 1, 1, &ok);
    int max_vote_threshold = QInputDialog::getInt(parent, "max_vote_threshold", "FIX ME:vote threshold (to be included for consensusing)", 0, 1, 1, 1, &ok);
    int cluster_distance_threshold = QInputDialog::getInt(parent, "FIX ME:cluster distance threshold", "5~10", 0, 1, 1, 1, &ok);

	if (!ok)
		return 0;
	
    QList<NeuronSWC> merge_result;



	QString fileSaveName;
	QString defaultSaveName = fileOpenName + "_consensus.swc";
	fileSaveName = QFileDialog::getSaveFileName(0, QObject::tr("Save merged neuron to file:"),
			defaultSaveName,
			QObject::tr("Supported file (*.swc)"
				";;Neuron structure	(*.swc)"
				));


    QString SelectedNeuronsAnoFileName = fileSaveName+"_SelectedNeurons.ano";
    //remove_outliers(nt_list, SelectedNeuronsAnoFileName);
    //if (!consensus_skeleton_votemap(nt_list, merge_result, max_vote_threshold, cluster_distance_threshold,  callback))
    if (!consensus_skeleton_match_center(nt_list, merge_result, max_vote_threshold,
                                         cluster_distance_threshold,  callback))
    {
        v3d_msg("Error in consensus swc!");
        return -1;
    }
    if (!export_listNeuron_2swc(merge_result,qPrintable(fileSaveName)))
	{
		v3d_msg("Unable to save file");
		return -1;
	}
    v3d_msg("Consensus swc is saved into: " + fileSaveName + "\n");
	return 1;
}

bool vote_map_func(const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 &callback)
{

    if(input.size()==0 || output.size() != 1) return false;
    char * paras = NULL;

    //parsing input
    vector<char *> * inlist =  (vector<char*> *)(input.at(0).p);
    if (inlist->size()==0)
    {
        cerr<<"You must specify input linker or swc files"<<endl;
        return false;
    }

    //parsing output
    vector<char *> * outlist = (vector<char*> *)(output.at(0).p);
    if (outlist->size()>1)
    {
        cerr << "You cannot specify more than 1 output files"<<endl;
        return false;
    }

    vector<NeuronTree> nt_list;
    QStringList nameList;
    QString qs_linker;
    char * dfile_result = NULL;
    V3DLONG neuronNum = 0;

    for (int i=0;i<inlist->size();i++)
    {
        qs_linker = QString(inlist->at(i));
        if (qs_linker.toUpper().endsWith(".ANO"))
        {
            cout<<"(0). reading a linker file."<<endl;
            P_ObjectFileType linker_object;
            if (!loadAnoFile(qs_linker,linker_object))
            {
                fprintf(stderr,"Error in reading the linker file.\n");
                return 1;
            }
            nameList = linker_object.swc_file_list;
            neuronNum += nameList.size();
            for (V3DLONG i=0;i<neuronNum;i++)
            {
                NeuronTree tmp = readSWC_file(nameList.at(i));
                nt_list.push_back(tmp);
            }
        }
        else if (qs_linker.toUpper().endsWith(".SWC"))
        {
            cout<<"(0). reading an swc file"<<endl;
            NeuronTree tmp = readSWC_file(qs_linker);
            nt_list.push_back(tmp);
            neuronNum++;
            if (outlist->size()==0)
            {
                cerr<<"You must specify outfile name if you input a list of swcs"<<endl;
                return false;
            }
        }
    }

    QString outfileName;
    if (outlist->size()==0)
        outfileName = qs_linker + "_vote_map.v3draw";
    else
        outfileName = QString(outlist->at(0));

    int dilation_radius = 0;
    if (!vote_map(nt_list,dilation_radius,outfileName, callback))
    {
        cerr<<"error in consensus_skeleton"<<endl;
        return false;
    }
  return true;
}

bool consensus_swc_func(const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 &callback)
{
    if(input.size()==0 || output.size() != 1) return false;
	char * paras = NULL;

	//parsing input
	vector<char *> * inlist =  (vector<char*> *)(input.at(0).p);
	if (inlist->size()==0)
	{
		cerr<<"You must specify input linker or swc files"<<endl;
		return false;
	}

	//parsing output
	vector<char *> * outlist = (vector<char*> *)(output.at(0).p);
	if (outlist->size()>1)
	{
        cerr << "You cannot specify more than 1 output files"<<endl;
		return false;
	}

	//parsing parameters

    int max_vote_threshold = 4;
    int cluster_distance_threshold = 10; // ignore nodes that far away for clustering
	if (input.size()==2)
	{
		vector<char*> * paras = (vector<char*> *)(input.at(1).p);
        if (paras->size() >= 1)
		{
            max_vote_threshold = atoi(paras->at(0));
            cout<<"max_vote_threshold = "<<max_vote_threshold<<endl;
            if (paras->size() == 2){
                cluster_distance_threshold =  atoi(paras->at(1));
                cout<<"clustering distance threshold = "<<cluster_distance_threshold<<endl;
            }
		}
		else
		{
			cerr<<"Too many parameters"<<endl;
			return false;
		}
	}

	vector<NeuronTree> nt_list;
	QStringList nameList;
	QString qs_linker;
	char * dfile_result = NULL;
	V3DLONG neuronNum = 0;

	for (int i=0;i<inlist->size();i++)
	{
		qs_linker = QString(inlist->at(i));
		if (qs_linker.toUpper().endsWith(".ANO"))
		{
			cout<<"(0). reading a linker file."<<endl;
			P_ObjectFileType linker_object;
            if (!loadAnoFile(qs_linker,linker_object))
			{
				fprintf(stderr,"Error in reading the linker file.\n");
				return 1;
			}
			nameList = linker_object.swc_file_list;
			neuronNum += nameList.size();
			for (V3DLONG i=0;i<neuronNum;i++)
			{
				NeuronTree tmp = readSWC_file(nameList.at(i));
				nt_list.push_back(tmp);
			}
		}
		else if (qs_linker.toUpper().endsWith(".SWC"))
		{
			cout<<"(0). reading an swc file"<<endl;
			NeuronTree tmp = readSWC_file(qs_linker);
			nt_list.push_back(tmp);
			neuronNum++;
			if (outlist->size()==0)
			{
				cerr<<"You must specify outfile name if you input a list of swcs"<<endl;
				return false;
			}
		}
	}

	QString outfileName;
	if (outlist->size()==0)
		outfileName = qs_linker + "_consensus.swc";
	else
		outfileName = QString(outlist->at(0));

	QList<NeuronSWC> merge_result;

    QString SelectedNeuronsAnoFileName = outfileName+"_SelectedNeurons.ano";
    // remove_outliers(nt_list, SelectedNeuronsAnoFileName);
    //if (!consensus_skeleton_vote_map(nt_list, merge_result, max_vote_threshold,cluster_distance_threshold, callback))
    if (!consensus_skeleton_match_center(nt_list, merge_result, max_vote_threshold,cluster_distance_threshold, callback))
	{
		cerr<<"error in consensus_skeleton"<<endl;
		return false;
	}

	export_listNeuron_2swc(merge_result,qPrintable(outfileName));
	printf("%s has been generated successfully\n",qPrintable(outfileName));

	return true;
}



bool dark_pruning_func(const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 &callback)
{

        if(input.size() <2 )
        {
            cerr<<"You must specify input eswc file and the corresponding image file."<<endl;
            return false;
        }

        //parsing input
        vector<char *> * inlist =  (vector<char*> *)(input.at(0).p);
        if (inlist->size()<2)
        {
            cerr<<"You must specify input eswc file and the corresponding image file."<<endl;
            return false;
        }

        //parsing output
        vector<char *> * outlist = (vector<char*> *)(output.at(0).p);
        if (outlist->size() > 1)
        {
            cerr << "You can only specify one output file"<<endl;
            return false;
        }

        //parsing parameters
        int visible_thre =0;
        vector<char*> * paras = (vector<char*> *)(input.at(1).p);
        if (paras->size() == 1)
        {
            visible_thre = atoi(paras->at(0));
            cout<<"The visible threshold for dark pruning is: "<<visible_thre<<endl;

        }
        else
        {
            cerr<<"One ( and only one) parameter is required."<<endl;
            return false;
        }



        QString input_swc_fn = QString(inlist->at(0));
        NeuronTree input_nt = readSWC_file(input_swc_fn);


         Image4DSimple * p4dImage = callback.loadImage( inlist->at(1) );
        if (!p4dImage || !p4dImage->valid())
            return false;



        QString outfileName;
        if (outlist->size()==0)
            outfileName = input_swc_fn + "_dark_pruned.swc";
        else
            outfileName = QString(outlist->at(0));

        QList<NeuronSWC> result_swc;
        if (!dark_pruning (input_nt, result_swc,p4dImage,visible_thre))
        {
            cerr<<"Error in dark pruning."<<endl;
            return false;
        }

        export_listNeuron_2swc(result_swc,qPrintable(outfileName));
        printf("%s has been generated successfully\n",qPrintable(outfileName));

        return true;


}





bool median_swc_func(const V3DPluginArgList & input, V3DPluginArgList & output)
{
    if(input.size()==0) return false;

    //parsing input
    vector<char *> * inlist =  (vector<char*> *)(input.at(0).p);
    if (inlist->size() == 0)
    {
        cerr<<"You must specify input linker or swc files"<<endl;
        return false;
    }

    //parsing output
    vector<char *> * outlist = (vector<char*> *)(output.at(0).p);
    if (outlist->size()>1)
    {
        cerr<<"You cannot specify more than 1 output files"<<endl;
        return false;
    }


    V3DLONG neuronNum = 0;
    vector<NeuronTree> nt_list;
    QString qs_linker;
    QStringList nameList;
    for (int i=0;i<inlist->size();i++)
    {
        qs_linker = QString(inlist->at(i));
        if (qs_linker.toUpper().endsWith(".ANO"))
        {
            cout<<"(0). reading a linker file."<<endl;
            P_ObjectFileType linker_object;
            if (!loadAnoFile(qs_linker,linker_object))
            {
                fprintf(stderr,"Error in reading the linker file.\n");
                return 1;
            }
            nameList = linker_object.swc_file_list;
            neuronNum += nameList.size();
            for (V3DLONG i=0;i<neuronNum;i++)
            {
                NeuronTree tmp = readSWC_file(nameList.at(i));
                nt_list.push_back(tmp);
            }
        }
        else if (qs_linker.toUpper().endsWith(".SWC") || qs_linker.toUpper().endsWith(".ESWC"))
        {
            cout<<"(0). reading an swc file"<<endl;
            NeuronTree tmp = readSWC_file(qs_linker);
            nt_list.push_back(tmp);
            neuronNum++;
        }
    }



    QString outfileName;
    if (outlist->size()==0)
        outfileName = qs_linker+ "_sum_dist.csv";
    else
        outfileName = QString(outlist->at(0));


    cout << "There are "<<nt_list.size() <<" input neurons."<<endl;

    //QString out_ano_file_name = outfileName + ".SelectedNeurons.ano";
    //remove_outliers(nt_list,out_ano_file_name);

    int idx = median_swc(nt_list,outfileName);
    if (idx <0){
        cerr << "error in median_swc()" << endl;
        return false;
    }
    QString fn = nt_list[idx].file;
    cout<<"Median swc is neuron " << idx <<" :" <<fn.toStdString().c_str()<< endl;
    return true;
}





bool average_node_position_func(const V3DPluginArgList & input, V3DPluginArgList & output)
{
    //parsing input
    vector<char *> * inlist =  (vector<char*> *)(input.at(0).p);
    cout<<"\n\n  inlist.size = "<<inlist->size()<<endl;
    if ( inlist->size() <2 )
    {
        cerr<<"You must specify inputs: median swc file and the linker file"<<endl;
        return false;
    }

    //parsing output
    vector<char *> * outlist = (vector<char*> *)(output.at(0).p);
    if (outlist->size()>1)
    {
        cerr<<"You cannot specify more than 1 output files"<<endl;
        return false;
    }

    //parsing parameters
    V3DLONG distance_threshold= 0;
    vector<char*> * paras = (vector<char*> *)(input.at(1).p);
    if (paras->size()==1)
    {
        distance_threshold = atoi(paras->at(0));
        cout<<"distance_threshold = "<<distance_threshold<<endl;
    }
    else
    {
        cerr<<"Too many parameters"<<endl;
        return false;
    }

    vector<NeuronTree> nt_list;
    QStringList nameList;
    QString qs_linker;
    char * dfile_result = NULL;
    V3DLONG neuronNum = 0;
    NeuronTree median_neuron;

    qs_linker = QString(inlist->at(0));
    if (qs_linker.toUpper().endsWith(".SWC"))
    {
        cout<<"(0). reading an swc file"<<endl;
        median_neuron= readSWC_file(qs_linker);
    }
    else{
        cout <<" The first input should be a swc file ( median neuron)." <<endl;
        return false;
    }

    qs_linker = QString(inlist->at(1));
    if (qs_linker.toUpper().endsWith(".ANO"))
    {
        cout<<"(0). reading a linker file."<<endl;
        P_ObjectFileType linker_object;
        if (!loadAnoFile(qs_linker,linker_object))
        {
            fprintf(stderr,"Error in reading the linker file.\n");
            return 1;
        }
        nameList = linker_object.swc_file_list;
        neuronNum += nameList.size();
        for (V3DLONG i=0;i<neuronNum;i++)
        {
            NeuronTree tmp = readSWC_file(nameList.at(i));
            nt_list.push_back(tmp);
        }
    }
    else
    {
        cout <<" The second input should be a ANO file ( the group of neurons)." <<endl;
    }

    QString outfileName;
    if (outlist->size()==0)
        outfileName = qs_linker + "_median_adjusted.swc";
    else
        outfileName = QString(outlist->at(0));

    NeuronTree median_adjusted = average_node_position(median_neuron, nt_list, distance_threshold );
    if (median_adjusted.listNeuron.size() == 0 ){

        cerr<<"error in average_node_position()"<<endl;
        return false;
    }

    writeSWC_file(outfileName, median_adjusted);
    printf("\t %s has been generated successfully\n",qPrintable(outfileName));

    return true;
}

int average_node_position_menu(V3DPluginCallback2 &callback, QWidget *parent)
{
    QString fileOpenName1;
    fileOpenName1 = QFileDialog::getOpenFileName(0, QObject::tr("Open the median SWC File"),
            "",
            QObject::tr("Supported file (*.swc)"
                ));
    if(fileOpenName1.isEmpty())
        return -1;
    NeuronTree median_neuron = readSWC_file(fileOpenName1);

     QString fileOpenName2 = QFileDialog::getOpenFileName(0, QObject::tr("Open the ano File that contains all input SWCs"),
            "",
            QObject::tr("Supported file (*.ano)"
                ";;Neuron structure	(*.ano)"
                ));
    if(fileOpenName2.isEmpty())
        return -1;

    P_ObjectFileType linker_object;
    if (!loadAnoFile(fileOpenName2,linker_object))
    {
        fprintf(stderr,"Error in reading the linker file.\n");
        return -1;
    }

    QStringList nameList = linker_object.swc_file_list;
    V3DLONG neuronNum = nameList.size();
    V3DLONG avg_node_num = 0;
    V3DLONG max_node_num = -1;
    vector<NeuronTree> nt_list;

    for (V3DLONG i=0;i<neuronNum;i++)
    {
        NeuronTree tmp = readSWC_file(nameList.at(i));
        nt_list.push_back(tmp);
        avg_node_num += tmp.listNeuron.size();
        if (tmp.listNeuron.size() > max_node_num)
            max_node_num = tmp.listNeuron.size();
    }
    avg_node_num /= neuronNum;

    bool ok;
    double distance_threshold= QInputDialog::getDouble(parent, "distance threshold",
              "Please specify the maximum distance allowed to search for maching nodes in all input neurons: ",8, 1, 100, 1, &ok );

    if (!ok)
        return 0;

    NeuronTree median_adjusted =  average_node_position(median_neuron, nt_list, distance_threshold);

    QString fileSaveName;
    QString defaultSaveName = fileOpenName1 + "_adjusted.swc";
    fileSaveName = QFileDialog::getSaveFileName(0, QObject::tr("Save adjusted median neuron to file:"),
            defaultSaveName,
            QObject::tr("Supported file (*.swc)"
                ";;Neuron structure	(*.swc)"
                ));
    if (!writeSWC_file(qPrintable(fileSaveName),median_adjusted))
    {
        v3d_msg("Unable to save file");
        return -1;
    }

    return 1;
}


int median_swc_menu(V3DPluginCallback2 &callback, QWidget *parent)
{
    QString fileOpenName;
    fileOpenName = QFileDialog::getOpenFileName(0, QObject::tr("Open the ANO File that contains all input SWC files"),
            "",
            QObject::tr("Supported file (*.ano)"
                ";;Neuron structure	(*.ano)"
                ));
    if(fileOpenName.isEmpty())
        return -1;

    P_ObjectFileType linker_object;
    if (!loadAnoFile(fileOpenName,linker_object))
    {
        fprintf(stderr,"Error in reading the linker file.\n");
        return -1;
    }

    QStringList nameList = linker_object.swc_file_list;
    V3DLONG neuronNum = nameList.size();
    vector<NeuronTree> nt_list;

    for (V3DLONG i=0;i<neuronNum;i++)
    {
        NeuronTree tmp = readSWC_file(nameList.at(i));
        nt_list.push_back(tmp);
    }

    v3d_msg( QString("It may take a while to compute pair-wise neuron distances. Please watch the output console for "
                     "computed metrics. Detailed distance reports are saved in "+fileOpenName+"_sum_dist.csv." ));


    QString outfileName = fileOpenName + "_sum_dist.csv";

    cout << "There are "<<nt_list.size() <<" input neurons."<<endl;

    QString out_ano_file_name = fileOpenName + ".SelectedNeurons.ano";
    //remove_outliers(nt_list,out_ano_file_name);

    int id = median_swc(nt_list,outfileName);

    if (id >=0 )
    {
        v3d_msg( QString("The median swc is : %1").arg( nameList.at(id)) );
    }

    /*
    QString fileSaveName;
    QString defaultSaveName = fileOpenName + "_consensus.swc";
    fileSaveName = QFileDialog::getSaveFileName(0, QObject::tr("Save merged neuron to file:"),
            defaultSaveName,
            QObject::tr("Supported file (*.swc)"
                ";;Neuron structure	(*.swc)"
                ));
    if ( !writeSWC_file(fiveSaveName, median_neuron))
    {
        v3d_msg("Unable to save file");
        return -1;
    }

    */
    return 1;
}

void printHelp()
{
    cout<<"\nConsensus Skeleton: This plugin has the following three functions:"<<endl;
    cout<<"\n  1) Pick the median neuron tree from a group of input neuron tress."<<endl;
    cout<<"\nUsage: v3d -x consensus_swc -f median_swc -i <input ANO linker file> [-o <output csv file>] "<<endl;
    cout<<"Parameters:"<<endl;
    cout<<"\t-f <function_name>:  median_swc"<<endl;
    cout<<"\t-i <input_file(s)>:  an input linker file (.ano) or multiple swc files"<<endl;
    cout<<"\t-o <output_csv_file>: print out the pair-wise distances for each input neuron to all other neurons."<<endl;
    cout<< " The index number of the median swc in the ano file will be reported in standard output. "<<endl;
    cout<<"Example: v3d -x consensus_swc -f median_swc -i mylinker.ano[myfolder/*.swc]  -o distances.csv \n"<<endl;

    cout<<"\n  2) Adjust input neuron node locations by averaging over the matching nodes from the input group of neurons tree."<<endl;
    cout<<"\nUsage: v3d -x consensus_swc -f average_node_position -i <median swc> <linker ANO file> -o <output_file> -p <distance_threshold>"<<endl;
    cout<<"Parameters:"<<endl;
    cout<<"\t-f <function_name>:  average_node_position"<<endl;
    cout<<"\t-i    <median swc>:  input median swc file (generated from median_swc function)"<<endl;
    cout<<"\t <linker ANO file>:  input linker file (.ano)"<<endl;
    cout<<"\t -p <distance_threshold>:  nodes that have distances larger than this threshold will "
          "not be considered matching for averaging." <<endl;
    cout<<"\t -o <output_file>:  output file name." <<endl;
    cout<<"Example: v3d -x consensus_swc -f average_node_position -i median.swc mylinker.ano -p 8 -o median_adjusted.swc\n"<<endl;

    cout<<"\n  3) Generate a consensus neuron skeleton (swc file) from a group of neurons ( radii are ignored)."<<endl;
    cout<<"\nUsage: v3d -x consensus_swc -f consensus_swc -i <input_folder or ano file> -o <output_file> "<<endl;
    cout<<"Parameters:"<<endl;
    cout<<"\t-f <function_name>:  consensus_swc"<<endl;
    cout<<"\t-i <input>:  input linker file (.ano) or folder path"<<endl;
    cout<<"\t-p <max_vote_threshold> <clustering_distance_threshold>:  max_vote_threshold: by default votes bigger than 1/3 of valid inputs will be "<<endl;
    cout<<"\t                                 included for consensing, this max_vote_threshold is setting the upper bound such voting threshold." <<endl;
    cout<<"\t-o <output_file>:  output file name. If -i is followd by a linker file name, this parameter can be omitted"<<endl;
    cout<<"\t                   default result will be generated under the same directory of the ref linkerfile and has a name of 'linkerFileName_consensus.swc'"<<endl;
    cout<<"Example: v3d -x consensus_swc -f consensus_swc -i mylinker.ano -o consensus.swc -p 6 5 \n"<<endl;
    cout<<"Example: v3d -x consensus_swc -f consensus_swc -i myfolder/*.swc -o consensus.swc -p 6 5 \n"<<endl;


    cout<<"\n  4) Generate a vote map volume (aggregated mask images) from multiple neurons ( radii are considered)."<<endl;
    cout<<"\nUsage: v3d -x consensus_swc -f vote_map -i <input> -o <output_image_file> "<<endl;
    cout<<"Parameters:"<<endl;
    cout<<"\t-f <function_name>:  vote_map"<<endl;
    cout<<"\t-i <input>:  input linker file (.ano) or folder path"<<endl;
    cout<<"\t-o <output_image_file>:  output image file name."<<endl;
    cout<<"Example: v3d -x consensus_swc -f vote_map -i mylinker.ano -o vote_map.v3draw\n"<<endl;

    cout<<"\n  5) Prune consensus based on original image" <<endl;
    cout<<"\nUsage: v3d -x consensus_swc -f dark_pruning -i <input_consensus_file> <input_image_file> -o <output_swc_file> "<<endl;
    cout<<"Parameters:"<<endl;
    cout<<"\t-f <function_name>:  dark_pruning"<<endl;
    cout<<"\t-i <input_eswc>  <input_image>:  input_consensus_eswc  input_image"<<endl;
    cout<<"\t-o <output_image_file>:  output image file name."<<endl;
    cout<<"\t-p <visible_threshold>:  visible threshold for dark pruning."<<endl;
    cout<<"Example: v3d -x consensus_swc -f dark_pruning -i input_consensus_file  input_image.v3dpbd -o pruned.swc -p 40\n"<<endl;
}


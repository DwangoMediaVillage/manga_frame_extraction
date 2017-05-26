//
//  main.cpp
//  MangaFrameExtraction
//
//  Created by 山田　祐雅　 on 2015/10/19.
//  Copyright (c) 2015年 山田　祐雅　. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <dirent.h>

#include <boost/regex.hpp>
#include <boost/program_options.hpp>

#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "fsp.h"

using namespace std;
using namespace boost::program_options;
using namespace cv;

int MAX_WIDTH = 1000;
int MAX_HEIGHT = 1000;

void execute(string input_dir, string output_dir, string filename)
{
    IplImage* input = cvLoadImage((input_dir + filename).c_str(), CV_LOAD_IMAGE_GRAYSCALE);

    // 幅・高さが最大長に収まるようにリサイズ。幅・高さがともに最大長より小さい場合は、ギリギリまで拡大
    if (input->width >= MAX_WIDTH || input->height >= MAX_HEIGHT) {
        int width = MAX_WIDTH;
        int height = (int)(input->height * ((double)MAX_WIDTH / input->width));
        if (height > MAX_HEIGHT) {
            width = (int)(width * ((double)MAX_HEIGHT / height));
            height = MAX_HEIGHT;
        }

        cout << "rescale" << " ";
        cout << (float) width / input->width << " ";
        cout << (float) height / input->height << " ";
        cout << endl;

        IplImage* scale_input = cvCreateImage(cvSize(width, height), input->depth, input->nChannels);
        cvResize(input, scale_input, CV_INTER_CUBIC);
        input = cvCloneImage(scale_input);
        cvReleaseImage(&scale_input);
    }

    PixPoint p = {0, 0};
    FrameSeparation* fs = new FrameSeparation(*input, filename, output_dir, input->width * input->height, p);
    cvReleaseImage(&input);
    separate_count = 0;
    delete fs;
}

int main(int argc, char* argv[])
{
    options_description mOption("Mode Option");
    options_description ioOption("I/O Option");
    options_description fileOption("File Convert Option");
    options_description helpOption("Help");

    mOption.add_options()("RecursiveMode,r", "recursive directory scan mode")("SingleMode,s", "single file scan mode");
    ioOption.add_options()("InputDir,d", value<string>(), "input directory (full path)")("InputFile,f", value<string>(), "input file (full path)")("OutputDir,o", value<string>(), "output directory (full path)");
    fileOption.add_options()("MaxWidth,w", value<int>()->default_value(1000), "max width")("MaxHeight,h", value<int>()->default_value(1000), "max height");
    helpOption.add_options()("Help,h", "view help");

    mOption.add(ioOption).add(fileOption).add(helpOption);
    variables_map values;

    try {
        store(parse_command_line(argc, argv, mOption), values);
        notify(values);
        MAX_WIDTH = values["MaxWidth"].as<int>();
        MAX_HEIGHT = values["MaxHeight"].as<int>();

        if (values.count("RecursiveMode") && values.count("SingleMode")) {
            cout << mOption << endl;
        }

        // recursive mode
        else if (values.count("RecursiveMode") && values.count("InputDir") && values.count("OutputDir")) {
            string input_dir = values["InputDir"].as<string>();
            string output_dir = values["OutputDir"].as<string>();

            DIR* dp;
            dirent* entry;
            boost::regex re(".*(.jpg).*");
            dp = opendir(input_dir.c_str());
            if (dp == NULL)
                exit(1);
            do {
                entry = readdir(dp);
                if (entry != NULL) {
                    if (regex_match(entry->d_name, re)) {
                        //#ifdef DEBUG
                        std::cout << entry->d_name << std::endl;
                        //#endif
                        //                execute(dir, output_dir, entry->d_name);
                        execute(input_dir, output_dir, entry->d_name);
                        //                exit(0);
                    }
                }
            } while (entry != NULL);
        }

        // single file mode
        else if (values.count("SingleMode") && values.count("InputFile") && values.count("OutputDir")) {
            string input_file = values["InputFile"].as<string>();
            string output_dir = values["OutputDir"].as<string>();
            string filename, path;

            size_t pos;
#ifdef _WIN32
            pos = input_file.find_last_of("\\");
#else
            pos = input_file.find_last_of("/");
#endif
            if (pos != std::string::npos) {
                path.assign(input_file.begin(), input_file.begin() + pos + 1);
                filename.assign(input_file.begin() + pos + 1, input_file.end());
                cout << "input file path:" << path << ", filename:" << filename << endl;
                execute(path, output_dir, filename);
            }
            else {
                cout << mOption << endl;
            }
        }

        else {
            cout << mOption << endl;
        }
    }
    catch (exception& e) {
        cout << e.what() << endl;
    }

    cout << "end" << endl;
}

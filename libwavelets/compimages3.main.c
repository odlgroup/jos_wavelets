#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#define COLORS 3
#define XLENGTH 512
#define YLENGTH 512
#define SIZE (COLOR * XLENGTH * YLENGTH)
#define LEVELS 6
#define SKIP 1
#define COMPARESKIP 1
#define max(s, t) (s > t ? s : t)
#define min(s, t) (s < t ? s : t)
#define TRIG_UNSIGNEDCHAR
#include "bio_parameters.h"
#ifndef FLOAT
#define FLOAT float
#endif

int compimages();
int show_compresults();

int main(argc, arglst) int argc;
char* arglst[];
{
    OUTTYPE* outvector;
    INTYPE* invector;
    FLOAT* tempvector;
    DIFFTYPE* diffvector = NULL;
    int optioncount = 0;
    int inputerror = 0;
    int argind;
    int Headerlength1 = 0;
    int Headerlength2 = 0;
    FILE* image1, *image2;
    int filerror = 0;
    FILE* argumentfile;
    int k;
    int Colors = 1;
    int Boxes = 1;
    int Xlength = XLENGTH;
    int Ylength = YLENGTH;
    int Zlength = ZLENGTH;
    long time00, time0, time1;
    int ifnotSilent = 1;

    int nrread;
    double maxvalue = 0;
    double poserror = 0;
    double negerror = 0;
    double error2 = 0;
    double norm2 = 0;
    double error1 = 0;
    double norm1 = 0;
    double dsum = 0;
    double errorlog10 = 0;
    int nonexactnumber = 0;
    int size = 0;
    int Size = 0;
    float SNR;
    time00 = clock();
    time0 = clock();

    if (!((argumentfile = fopen("comp.arg", "r")) == NULL)) {
        Headerlength1 = _getw(argumentfile);
        Headerlength2 = _getw(argumentfile);
        Boxes = _getw(argumentfile);
        Colors = _getw(argumentfile);
        Xlength = _getw(argumentfile);
        Ylength = _getw(argumentfile);
        Zlength = _getw(argumentfile);
        ifnotSilent = _getw(argumentfile);
        fclose(argumentfile);
    }

    for (argind = 1; argind < argc; argind++) {
        if (strcmp(arglst[argind], "-Header1") == 0) {
            if ((argind < 1) || (argind == argc - 1))
                inputerror = 1;
            else
                Headerlength1 = atoi(arglst[argind + 1]);
            optioncount += 2;
        }
    }

    for (argind = 1; argind < argc; argind++) {
        if (strcmp(arglst[argind], "-Header2") == 0) {
            if ((argind < 1) || (argind == argc - 1))
                inputerror = 1;
            else
                Headerlength2 = atoi(arglst[argind + 1]);
            optioncount += 2;
        }
    }

    for (argind = 1; argind < argc; argind++) {
        if ((strcmp(arglst[argind], "-Color") && strcmp(arglst[argind], "-Grey")) == 0) {
            if ((argind < 1) || (argind > argc - 1))
                inputerror = 1;
            else {
                if (strcmp(arglst[argind], "-Color") == 0)
                    Colors = 3;
                else
                    Colors = 1;
            }
            optioncount += 1;
        }
    }

    for (argind = 1; argind < argc; argind++) {
        if (strcmp(arglst[argind], "-L") == 0) {
            if ((argind < 1) || (argind >= argc - 3))
                inputerror = 1;
            else {
                Xlength = atoi(arglst[argind + 1]);
                Ylength = atoi(arglst[argind + 2]);
                Zlength = atoi(arglst[argind + 3]);
            }
            optioncount += 4;
        }
    }

    for (argind = 1; argind < argc; argind++) {
        if (strcmp(arglst[argind], "-Boxes") == 0) {
            if ((argind < 1) || (argind == argc - 1))
                inputerror = 1;
            else
                Boxes = atoi(arglst[argind + 1]);
            optioncount += 2;
        }
    }

    for (argind = 1; argind < argc; argind++) {
        if ((strcmp(arglst[argind], "-Silent") && strcmp(arglst[argind], "-NotSilent")) == 0) {
            if ((argind < 1) || (argind > argc - 1))
                inputerror = 1;
            else {
                if (strcmp(arglst[argind], "-NotSilent") == 0)
                    ifnotSilent = 1;
                else
                    ifnotSilent = 0;
            }
            optioncount += 1;
        }
    }

    if (inputerror || ((argc > 3 + optioncount) || (argc < 1 + optioncount))) {
        printf("Error: Arguments:  3  + %d options  but  argc = %d.  inputerror=%d\n", optioncount, argc, inputerror);
        inputerror = 2;
    } else {
        if (argc < 3 + optioncount) inputerror = 1;
    }

    if (inputerror > 0) {
        printf("use:  compimages\n\t [<filename1><filename2>\n");
        printf("\t\t(without this option only parameters will be changed)]\n");
        printf("\t[-Header1 <header length>]\n");
        printf("\t[-Header2 <header length>]\n");
        printf("\t[-Color || -Grey ]\n");
        printf("\t[-Boxes <number of boxes>]\n");
        printf("\t[-L <xlength> <ylength> <zlength>]\n");
        printf("\t[ -NotSilent || -Silent] \n");
    }

    if (inputerror > 1) return 2;
    if (ifnotSilent) {
        printf("Comparing two file with %d boxes: (%d x %d x %d in %d colors):\n",
               Boxes, Xlength, Ylength, Zlength, Colors);
        printf("First %d bytes in <filename1> are skipped as header.\n", Headerlength1);
        printf("First %d bytes in <filename2> will be skipped as header.\n\n", Headerlength2);
    }
    argumentfile = fopen("comp.arg", "w");
    _putw(Headerlength1, argumentfile);
    _putw(Headerlength2, argumentfile);
    _putw(Boxes, argumentfile);
    _putw(Colors, argumentfile);
    _putw(Xlength, argumentfile);
    _putw(Ylength, argumentfile);
    _putw(Zlength, argumentfile);
    _putw(ifnotSilent, argumentfile);
    fclose(argumentfile);

    if (inputerror > 0) return 1;

    size = Boxes * Colors * Ylength * Xlength * Zlength;
    invector = (INTYPE*)malloc((unsigned int)size * sizeof(INTYPE));
    outvector = (OUTTYPE*)malloc((unsigned int)size * sizeof(OUTTYPE));
    /*diffvector=(DIFFTYPE*)malloc((unsigned int)size*sizeof(DIFFTYPE));*/
    tempvector = (FLOAT*)malloc((unsigned int)((size)) * sizeof(FLOAT));
    if ((image1 = fopen(arglst[1], "r")) == NULL) {
        printf("Cannot open file %s\n", arglst[1]);
        filerror = 1;
    }
    if ((image2 = fopen(arglst[2], "r")) == NULL) {
        printf("Cannot open file %s\n", arglst[2]);
        filerror = 1;
    }
    if (filerror) return 1;
    if (ifnotSilent)
        printf("Loading image1\n");
    time0 = clock();

    if (0 && Headerlength1) {
        if (ifnotSilent)
            printf("Reading %d first bytes in %s,as header\n", Headerlength1, arglst[1]);
        for (k = 0; k < Headerlength1; k++)
            fgetc(image1);
    }
    nrread = fread(invector, sizeof(INTYPE), size, image1);
    if (nrread != size) {
        printf("Error reading invector.\nOnly %d types of %d are read\n",
               nrread, size);
        filerror = 1;
    }
    fclose(image1);

    printf("sizeof(INTYPE)=%lu, sizeof(OUTTTYPE)%lu\n", sizeof(INTYPE), sizeof(OUTTYPE));
    time1 = clock();
    if (ifnotSilent) {
        printf("time used: %f\n", ((float)(time1 - time0) / CLOCKS_PER_SEC));
        printf("Loading image2:\n");
    }
    time0 = clock();
    image2 = fopen(arglst[2], "r");
    if (0 && Headerlength2) {
        if (ifnotSilent)
            printf("Reading %d first bytes in %s as header\n", Headerlength2, arglst[2]);
        for (k = 0; k < Headerlength2; k++)
            fgetc(image2);
    }
    nrread = fread(outvector, sizeof(INTYPE), size, image2);
    if (nrread != size) {
        printf("Error reading outvector.\nOnly %d types of %d are read\n",
               nrread, size);
        filerror = 1;
    }
    fclose(image2);

    time1 = clock();
    if (ifnotSilent) {
        printf("time used: %f\n", ((float)(time1 - time0) / CLOCKS_PER_SEC));
        printf("\nComparing with original:\n");
    }

    compimages(invector, outvector, diffvector, size,
               &maxvalue, &poserror, &negerror, &error2, &norm2, &error1, &norm1,
               &errorlog10, &nonexactnumber, &Size, &dsum);

    free(invector);
    free(outvector);
    free(tempvector);

    show_compresults(maxvalue, poserror, negerror, error2, norm2, error1, norm1,
                     errorlog10, nonexactnumber, size, &SNR, sizeof(INTYPE), dsum);
    /*
if(ifnotSilent)
  printf("\nWriting differenses to file: %s\n", "diff_fileH");
diffile=fopen("diff_fileH","w");
if(Colors==3){ FPRINTF_COLORHEADER(diffile,Xlength,Ylength);}
if(Colors==1){FPRINTF_GREYHEADER(diffile,Xlength,Ylength);}

fwrite(diffvector,sizeof(DIFFTYPE),size,diffile);
fclose(diffile);
*/
    return 0;
}

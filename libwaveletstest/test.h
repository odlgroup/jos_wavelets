# Copyright 2015-2016 Jan-Olov Stromberg <jostromb@kth.se>
#
# jos_wavelets is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# jos_wavelets is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with jos_wavelets.  If not, see <http://www.gnu.org/licenses/>.



int compimages(
    FLOAT*,  //vector1
    FLOAT*,  //vector2
    FLOAT*,  //diffvector
    int,     //Size
    double*, //ptrmaxvalue
    double*, //ptrposerror
    double*, //ptrnegerror
    double*, //ptrerror2
    double*, //ptrnorm2
    double*, //ptrerror1
    double*, //ptrnorm1
    double*, //ptrerrorlog10
    int*,    //ptrnonexactnumber
    int*,    //ptrsize
    double*  //ptrsum
    );

int show_compresults(double, // maxvalue
                     double, // poserror
                     double, // negerror
                     double, // error2
                     double, // norm2
                     double, // error1
                     double, // norm1
                     double, // errorlog10
                     int,    // nonexactnumber
                     int,    // size
                     FLOAT*, // *SNR_ptr
                     int,    // typesize
                     double  // sum
                     );

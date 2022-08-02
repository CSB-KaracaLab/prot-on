/*******************************************************************************************************************************
This file is a part of the EvoDesign physical Energy Function (EvoEF)

Copyright (c) 2019 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

#include "ErrorHandling.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

BOOL FAILED(int errorCode){
  if(errorCode == Success || errorCode == Warning){
    return FALSE;
  }
  return TRUE;
}

int TraceError(char* userMsg, int errorCode){
  char errMsg[MAX_LENGTH_ERR_MSG+1];

  if(errorCode == Success){
    return errorCode;
  }
  else if (errorCode == Warning){
    strcpy(errMsg, "Warning message,");
    printf("--------------------------------------------\n");
    printf("Warning %d : %s %s.\n", errorCode, errMsg, userMsg);
    return errorCode;
  }

  switch(errorCode){
    case IOError:
      strcpy(errMsg, "File cannot be opened,"); break;
    case FormatError:
      strcpy(errMsg, "File format is wrong,"); break;
    case IndexError:
      strcpy(errMsg, "Index not in range,"); break;
    case ValueError:
      strcpy(errMsg, "Invalid value used,"); break;
    case ZeroDivisonError:
      strcpy(errMsg, "Denominator is near zero,"); break;
    case DataNotExistError:
      strcpy(errMsg, "Data not exist in records,"); break;
    case NameError:
      strcpy(errMsg, "Parameter name is wrong,"); break;
    default:
      strcpy(errMsg, "Undefined error type,");
  }

  printf("--------------------------------------------\n");
  printf("Error %d : %s %s.\n", errorCode, errMsg, userMsg);
  return errorCode;
}

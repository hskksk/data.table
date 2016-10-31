#include "data.table.h"
#include <errno.h>
#include <unistd.h>  // for access()
#include <fcntl.h>
#include <time.h>
#ifdef WIN32
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#define WRITE _write
#define CLOSE _close
#else
#define WRITE write
#define CLOSE close
#endif

static inline int maxStrLen(SEXP x, int na_len) {
  // The max nchar of any string in a column or factor level
  int max=na_len, nrow=length(x);
  for (int i=0; i<nrow; i++) {
    int l = LENGTH(STRING_ELT(x,i));
    if (l>max) max=l;
  }
  // TODO if(quote) count them here (we hopped already to get LENGTH), add that quote count to max
  //      exactly and remove the *(1+quote) later
  //      if there are no quotes present in any string in a column save in quotePresent[ncol] lookup
  //      and switch to memcpy when it's known no quotes are present to be escaped for that column
  return max;
}

#define QUOTE_FIELD \
  *ch++ = QUOTE; \
  for (const char *ch2 = CHAR(str); *ch2 != '\0'; ch2++) { \
    if (*ch2 == QUOTE) *ch++ = ESCAPE_QUOTE; \
    if (qmethod_escape && *ch2 == ESCAPE) *ch++ = ESCAPE; \
    *ch++ = *ch2; \
  } \
  *ch++ = QUOTE

#define NUM_SF 15
#define DECIMAL_SEP '.'  // TODO allow other decimal separator e.g. ','

// Globals for this file only (written once to hold parameters passed from R level)
static int na_len;
static const char *na_str;

static inline void writeIntegerCore(int x, char **thisCh)
{
  char *ch = *thisCh;
  if (x == NA_INTEGER) {
    if (na_len) { memcpy(ch, na_str, na_len); ch += na_len; }
  } else if (x == 0) {
    *ch++ = '0';
  } else {
    if (x<0) { *ch++ = '-'; x=-x; }
    // avoid log() call for speed. write backwards then reverse when we know how long
    int width = 0;
    while (x>0) { *ch++ = '0'+x%10; x /= 10; width++; }
    for (int i=width/2; i>0; i--) {
      char tmp=*(ch-i);
      *(ch-i) = *(ch-width+i-1);
      *(ch-width+i-1) = tmp;
    }
  }
  *thisCh = ch;
}

static inline void writeInteger(int col, int x, char **thisCh)
{
  char *ch = *thisCh;
  if (x == NA_INTEGER) {
    //do nothing
  } else {
    if(col==0){
      writeIntegerCore(x, &ch);
    }else{
      writeIntegerCore(col, &ch);
      *ch++ = ':';
      writeIntegerCore(x, &ch);
    }
  }
  *thisCh = ch;
}

union {
  double d;
  unsigned long long ull;
} u;

extern Rboolean verbose; // set by writefile

static inline void writeNumeric(int col, double x, char **thisCh)
{
  // hand-rolled / specialized for speed
  // *thisCh is safely the output destination with enough space (ensured via calculating maxLineLen up front)
  // technique similar to base R (format.c:formatReal and printutils.c:EncodeReal0)
  // differences/tricks :
  //   i) no buffers. writes straight to the final file buffer passed to write()
  //  ii) no C libary calls such as sprintf() where the fmt string has to be interpretted over and over
  // iii) no need to return variables or flags.  Just writes.
  //  iv) shorter, easier to read and reason with. In one self contained place.
  char *ch = *thisCh;
  // write index no when x is not NA
  if (R_FINITE(x) || !ISNAN(x)) {
    if(col>0){
      writeIntegerCore(col, &ch);
      *ch++ = ':';
    }
  }
  if (!R_FINITE(x)) {
    if (ISNAN(x)) {
      //do nothing
    } else if (x>0) {
      *ch++ = 'I'; *ch++ = 'n'; *ch++ = 'f';
    } else {
      *ch++ = '-'; *ch++ = 'I'; *ch++ = 'n'; *ch++ = 'f';
    }
  } else if (x == 0.0) {
    *ch++ = '0';   // and we're done.  so much easier rather than passing back special cases
  } else {
    if (x < 0.0) { *ch++ = '-'; x = -x; }  // and we're done on sign, already written. no need to pass back sign
    
    /*
    int exp = (int)floorl(log10l((long double)x));
    unsigned long long l = (unsigned long long)((long double)x * powl(10.0L, NUM_SF-exp));
    */
    
    u.d = x;
    unsigned long long fraction = u.ull & ((1ULL<<52)-1);  // TODO const & 0x...
    int exponent = ((int)(u.ull>>52) & 0x7FF) - 1023;
    long double acc = 0;
    for (int i=0; i<52; i++) { 
      // important smallest first; i.e. 2^-52
      if ( fraction & (1ULL<<i) ) acc+=ldexpl(1.0L,i-52);
      // TODO lookup table and avoid branch
      // TODO skip batches of 0 in the fraction by skipping over trailing 0 bytes perhaps
    }
    long double y = ldexpl(1.0L+acc, exponent);
    
    //long double y = 1.0 + fraction / 4503599627370496.0L;  // 2^52 as long double
    //y = ldexpl(y, exponent-1023);
    // long double y = (long double)frexpl(x, &e2);
    // int nd = (int)(e2 * 0.30102999566);  // * log10(2)
    // y *= (powl(2.0, e2) / powl(10.0, nd));  // doesn't seem efficient or accurate
    
    // since y was 1.* above so surely we know
    int exp = (int)floorl(log10l(y));
    unsigned long long l = (unsigned long long)(y * powl(10.0L, NUM_SF-exp));
    
    if (verbose) Rprintf("\nTRACE: acc=%.20Le ; y=%.20Le ; l=%llu ; e=%d     ", acc, y, l, exp);
    
    /*
    if (y<1) { y *= 1e15; exp++; }
    else  
    y *= 1e15;
    unsigned long long l = (unsigned long long)y;
    int exp = 
    
    int m = log10(2^e2)
    y * powl(10.0, NUM_DF)
    
    int exp = (int)floor(log10(x));
    
    unsigned long long l = (unsigned long long)(y *
                                                 (long double)powl((long double)10,(long double)(NUM_SF-exp)));
    //exp+=nd;
    */
    
    // TODO?: use lookup table like base R? ................. ^^^^
    //        here in fwrite it might make a difference whereas in base R other very
    //        significant write.table inefficiency dominates.
    // long double needed for 1729.9 to ensure 1e-310 doesn't write as 0.
    // Peppered with long double as windows seems to need that
    // l now contains NUM_SF+1 digits.
    // ULL for compound accuracy. If double, the repeated base 10 ops below could compound errors
    if (l%10 >= 5) l+=10; // use the last digit to round
    l /= 10;
    if (l == 0) {
      if (*(ch-1)=='-') ch--;
      *ch++ = '0';
    } else {
      // Count trailing zeros and therefore s.f. present in l
      int trailZero = 0;
      while (l%10 == 0) { l /= 10; trailZero++; }
      int sf = NUM_SF - trailZero;
      if (sf==0) {sf=1; exp++;}  // e.g. l was 9999999[5-9] rounded to 10000000 which added 1 digit
      
      // l is now an unsigned long that doesn't start or end with 0
      // sf is the number of digits now in l
      // exp is e<exp> were l to be written with the decimal sep after the first digit
      int dr = sf-exp-1; // how many characters to print to the right of the decimal place
      int width=0;       // field width were it written decimal format. Used to decide whether to or not.
      int dl0=0;         // how many 0's to add to the left of the decimal place before starting l
      if (dr<=0) { dl0=-dr; dr=0; width=sf+dl0; }  // 1, 10, 100, 99000
      else {
        if (sf>dr) width=sf+1;                     // 1.234 and 123.4
        else { dl0=1; width=dr+1+dl0; }            // 0.1234, 0.0001234
      }
      // So:  3.1416 => l=31416, sf=5, exp=0     dr=4; dl0=0; width=6
      //      30460  => l=3046, sf=4, exp=4      dr=0; dl0=1; width=5
      //      0.0072 => l=72, sf=2, exp=-3       dr=4; dl0=1; width=6
      if (width <= sf + (sf>1) + 2 + (abs(exp)>99?3:2)) {
         //              ^^^^ to not include 1 char for dec in -7e-04 where sf==1
         //                      ^ 2 for 'e+'/'e-'
         // decimal format ...
         ch += width-1;
         if (dr) {
           while (dr && sf) { *ch--='0'+l%10; l/=10; dr--; sf--; }
           while (dr) { *ch--='0'; dr--; }
           *ch-- = DECIMAL_SEP;
         }
         while (dl0) { *ch--='0'; dl0--; }
         while (sf) { *ch--='0'+l%10; l/=10; sf--; }
         // ch is now 1 before the first char of the field so position it afterward again, and done
         ch += width+1;
      } else {
        // scientific ...
        ch += sf;  // sf-1 + 1 for dec
        for (int i=sf; i>1; i--) {
          *ch-- = '0' + l%10;   
          l /= 10;
        }
        if (sf == 1) ch--; else *ch-- = DECIMAL_SEP;
        *ch = '0' + l;
        ch += sf + (sf>1);
        *ch++ = 'e';  // lower case e to match base::write.csv
        if (exp < 0) { *ch++ = '-'; exp=-exp; }
        else { *ch++ = '+'; }  // to match base::write.csv
        if (exp < 100) {
          *ch++ = '0' + (exp / 10);
          *ch++ = '0' + (exp % 10);
        } else {
          *ch++ = '0' + (exp / 100);
          *ch++ = '0' + (exp / 10) % 10;
          *ch++ = '0' + (exp % 10);
        }
      }
    }
  }
  *thisCh = ch;
}

static int NROWS(SEXP x){
  if(isMatrix(x)){
    return nrows(x);
  }else{
    return length(x);
  }
}

static int NCOLS(SEXP x){
  if(isMatrix(x)){
    return ncols(x);
  }else{
    return 1;
  }
}

SEXP writefile_libsvm(SEXP list_of_matrices,
                      SEXP filenameArg,
                      SEXP col_sep_Arg,
                      SEXP row_sep_Arg,
                      SEXP na_Arg,
                      SEXP quoteArg,           // TRUE|FALSE
                      SEXP qmethod_escapeArg,  // TRUE|FALSE
                      SEXP append,             // TRUE|FALSE
                      SEXP col_names,          // TRUE|FALSE
                      SEXP verboseArg,
                      SEXP turboArg)
{
  if (!isNewList(list_of_matrices)) error("fwrite must be passed an object of type list, data.table or data.frame");
  RLEN nmat = length(list_of_matrices);
  if (nmat==0) error("fwrite must be passed a non-empty list");

  // count the number of columns
  int ncols = 0;
  for (int i=1; i<nmat; i++) {
    SEXP mat = VECTOR_ELT(list_of_matrices, i);
    ncols += NCOLS(mat);
  }

  RLEN nrows = NROWS(VECTOR_ELT(list_of_matrices, 0));
  for (int i=1; i<nmat; i++) {
    if (nrows != NROWS(VECTOR_ELT(list_of_matrices, i)))
      error("Column %d's length (%d) is not the same as column 1's length (%d)", i+1, NROWS(VECTOR_ELT(list_of_matrices, i)), nrows);
  }
#ifndef _OPENMP
  Rprintf("Your platform/environment has not detected OpenMP support. fwrite() will still work, but slower in single threaded mode.\n");
  // Rprintf rather than warning() because warning() would cause test.data.table() to error about the unexpected warnings
#endif 
  verbose = LOGICAL(verboseArg)[0];
  const Rboolean quote = LOGICAL(quoteArg)[0];
  const Rboolean turbo = LOGICAL(turboArg)[0];
  
  const char col_sep = *CHAR(STRING_ELT(col_sep_Arg, 0));  // DO NOT DO: allow multichar separator (bad idea)
  
  const char *row_sep = CHAR(STRING_ELT(row_sep_Arg, 0));
  const int row_sep_len = strlen(row_sep);  // someone somewhere might want a trailer on every line
  na_str = CHAR(STRING_ELT(na_Arg, 0));
  na_len = strlen(na_str);
  const char QUOTE = '"';
  const char ESCAPE = '\\';
  const Rboolean qmethod_escape = LOGICAL(qmethod_escapeArg)[0];
  const char ESCAPE_QUOTE = qmethod_escape ? ESCAPE : QUOTE;
  const char *filename = CHAR(STRING_ELT(filenameArg, 0));

  errno = 0;   // clear flag possibly set by previous errors
  int f;
  if (*filename=='\0') f=-1;  // file="" means write to standard output
  else { 
#ifdef WIN32
    f = _open(filename, _O_WRONLY | _O_BINARY | _O_CREAT | (LOGICAL(append)[0] ? _O_APPEND : _O_TRUNC), _S_IWRITE);
    // row_sep must be passed from R level as '\r\n' on Windows since write() only auto-converts \n to \r\n
    // in _O_TEXT mode. We use O_BINARY for full control and perhaps speed since O_TEXT must have to deep branch an if('\n')
#else
    f = open(filename, O_WRONLY | O_CREAT | (LOGICAL(append)[0] ? O_APPEND : O_TRUNC), 0644);
#endif
    if (f == -1) {
      if( access( filename, F_OK ) != -1 )
        error("File exists and failed to open for writing. Do you have write permission to it? Is this Windows and does another process such as Excel have it open? File: %s", filename);
      else 
        error("Unable to create new file for writing (it does not exist already). Do you have permission to write here and is there space on the disk? File: %s", filename); 
    }
  }
  int true_false;
  
  clock_t t0=clock();
  // i) prefetch levels of factor columns (if any) to save getAttrib on every field on every row of any factor column
  // ii) calculate certain upper bound of line length
  SEXP levels[ncols];  // on-stack vla
  int lineLenMax = 2 + ncols*(int)(log10(ncols)+1);  // initialize with eol max width of \r\n on windows + length of indices
  int sameType = TYPEOF(VECTOR_ELT(list_of_matrices, 0));
  for (int mat_i=0; mat_i<nmat; mat_i++) {
    SEXP mat = VECTOR_ELT(list_of_matrices, mat_i);
    int mat_ncols = NCOLS(mat);
    if (TYPEOF(mat) != sameType) sameType = 0;
    for (int mat_col_i=0; mat_col_i<mat_ncols; mat_col_i++) {
      switch(TYPEOF(mat)) {
      case LGLSXP:
        lineLenMax+=1*mat_ncols;  // width of FALSE
        break;
      case REALSXP:
        lineLenMax+=25*mat_ncols;   // +- 15digits dec e +- nnn = 22 + 3 safety = 25
        break;
      case INTSXP:
        if (isFactor(mat)) {
          //TODO: implement case STRSXP
          error("Column %d's type is '%s' - not yet implemented.", mat_col_i+1,type2char(TYPEOF(mat)) );
          levels[mat_col_i] = getAttrib(mat, R_LevelsSymbol);
          sameType = 0; // TODO: enable deep-switch-avoidance for all columns factor
          lineLenMax += maxStrLen(levels[mat_col_i], na_len)*(1+quote) + quote*2;
          //                                    ^^^^^^^^^^ in case every character in the field is a quote, each to be escaped (!)
        } else {
          levels[mat_col_i] = NULL;
          lineLenMax+=11*mat_ncols;   // 11 + sign
        }
        break;
      case STRSXP:
        //TODO: implement case STRSXP
        error("Column %d's type is '%s' - not yet implemented.", mat_col_i+1,type2char(TYPEOF(mat)) );
        lineLenMax += maxStrLen(mat, na_len)*(1+quote) + quote*2;
        break;
      default:
        error("Column %d's type is '%s' - not yet implemented.", mat_col_i+1,type2char(TYPEOF(mat)) );
      }
      lineLenMax++;  // column separator
    }
  }
  clock_t tlineLenMax=clock()-t0;
  if (verbose) Rprintf("Maximum line length is %d calculated in %.3fs\n", lineLenMax, 1.0*tlineLenMax/CLOCKS_PER_SEC);
  // TODO: could parallelize by column, but currently no need as insignificant time
  t0=clock();

  if (nrows == 0) {
    if (verbose) Rprintf("No data rows present (nrow==0)\n");
    if (f!=-1 && CLOSE(f)) error("Error closing file: %s", filename);
    return(R_NilValue);
  }

  // Decide buffer size on each core
  // Large enough to fit many lines (to reduce calls to write()). Small enough to fit in each core's cache.
  // If the lines turn out smaller, that's ok. we just won't use all the buffer in that case. But we must be
  // sure to allow for worst case; i.e. every row in the batch all being the maximum line length.
  int bufSize = 1*1024*1024;  // 1MB  TODO: experiment / fetch cache size
  if (lineLenMax > bufSize) bufSize = lineLenMax;
  const int rowsPerBatch = bufSize/lineLenMax;
  const int numBatches = (nrows-1)/rowsPerBatch + 1;
  if (verbose) Rprintf("Writing data rows in %d batches of %d rows (each buffer size %.3fMB, turbo=%d) ... ",
    numBatches, rowsPerBatch, 1.0*bufSize/(1024*1024), turbo);
  t0 = clock();
  
  int nth;
  #pragma omp parallel num_threads(getDTthreads())
  {
    char *buffer = malloc(bufSize);  // one buffer per thread
    // TODO Ask Norm how to error() safely ... if (buffer == NULL) error("Unable to allocate %dMB buffer", bufSize/(1024*1024)); 
    char *ch = buffer;
    #pragma omp single
    {
      nth = omp_get_num_threads();
    }    
    #pragma omp for ordered schedule(dynamic)
    for(RLEN start_row = 0; start_row < nrows; start_row += rowsPerBatch) { 
      int upp = start_row + rowsPerBatch;
      if (upp > nrows) upp = nrows;
      if (turbo && sameType == REALSXP) {
        // avoid deep switch. turbo switches on both sameType and specialized writeNumeric
        for (RLEN row_i = start_row; row_i < upp; row_i++) {
          int col_i = 0;
          for (int mat_i = 0; mat_i < ncols; mat_i++) {
            SEXP mat = VECTOR_ELT(list_of_matrices, mat_i);
            int mat_ncols = NCOLS(mat);
            for (int mat_col_i = 0; mat_col_i < mat_ncols; mat_col_i++) {
              int mat_pos_i = row_i+mat_col_i*nrows;
              writeNumeric(col_i, REAL(mat)[mat_pos_i], &ch);
              *ch++ = col_sep;
              col_i++;
            }
          }
          ch--;  // backup onto the last col_sep after the last column
          memcpy(ch, row_sep, row_sep_len);  // replace it with the newline.
          ch += row_sep_len;
        }
      } else if (turbo && sameType == INTSXP) {
        for (RLEN row_i = start_row; row_i < upp; row_i++) {
          int col_i = 0;
          for (int mat_i = 0; mat_i < ncols; mat_i++) {
            SEXP mat = VECTOR_ELT(list_of_matrices, mat_i);
            int mat_ncols = NCOLS(mat);
            for (int mat_col_i = 0; mat_col_i < mat_ncols; mat_col_i++) {
              int mat_pos_i = row_i+mat_col_i*nrows;
              writeInteger(col_i, INTEGER(mat)[mat_pos_i], &ch);
              *ch++ = col_sep;
              col_i++;
            }
          }
          ch--;  // backup onto the last col_sep after the last column
          memcpy(ch, row_sep, row_sep_len);  // replace it with the newline.
          ch += row_sep_len;
        }
      } else {
        for (RLEN row_i = start_row; row_i < upp; row_i++) {
          int col_i = 0;
          for (int mat_i = 0; mat_i < nmat; mat_i++) {
            SEXP mat = VECTOR_ELT(list_of_matrices, mat_i);
            int mat_ncols = NCOLS(mat);
            for (int mat_col_i = 0; mat_col_i < mat_ncols; mat_col_i++) {
              int mat_pos_i = row_i+mat_col_i*nrows;
              SEXP str;
              switch(TYPEOF(mat)) {
              case LGLSXP:
                true_false = LOGICAL(mat)[mat_pos_i];
                if (true_false == NA_LOGICAL) {
                  if(col_i==0){
                    memcpy(ch, "NaN", 1);
                    ch += 1;
                  }
                } else if (true_false) {
                  if(col_i==0){
                    memcpy(ch, "1", 1);
                    ch += 1;
                  }else{
                    ch += sprintf(ch, "%d:1", col_i);
                  }
                } else {
                  if(col_i==0){
                    memcpy(ch, "0", 1);
                    ch += 1;
                  }else{
                    ch += sprintf(ch, "%d:0", col_i);
                  }
                }
                break;
              case REALSXP:
                if (ISNA(REAL(mat)[mat_pos_i])) {
                  if(col_i==0){
                    memcpy(ch, "NaN", 1);
                    ch += 1;
                  }
                } else {
                  if (turbo) {
                    // if there are any problems with the hand rolled double writing, then turbo=FALSE reverts to standard library
                    writeNumeric(col_i, REAL(mat)[mat_pos_i], &ch);
                  } else {
                    //tt0 = clock();
                    if(col_i==0){
                      ch += sprintf(ch, "%.15G", REAL(mat)[mat_pos_i]);
                    }else{
                      ch += sprintf(ch, "%d:%.15G", col_i, REAL(mat)[mat_pos_i]);
                    }
                    //tNUM += clock()-tt0;
                  }
                }
                break;
              case INTSXP:
                if (INTEGER(mat)[mat_pos_i] == NA_INTEGER) {
                  // do nothing
                  if(col_i==0){
                    memcpy(ch, "NaN", 1);
                    ch += 1;
                  }
                } else if (levels[mat_col_i] != NULL) {   // isFactor(mat) == TRUE
                  // TODO: implement factor case
                } else {
                  if (turbo) {
                    writeInteger(col_i, INTEGER(mat)[mat_pos_i], &ch);
                  } else {
                    if(col_i==0){
                      ch += sprintf(ch, "%d", INTEGER(mat)[mat_pos_i]);
                    }else{
                      ch += sprintf(ch, "%d:%d", col_i, INTEGER(mat)[mat_pos_i]);
                    }
                  }
                }
                break;
              case STRSXP:
                // TODO: implement case STRSXP
                //str = STRING_ELT(mat, i);
                //if (str==NA_STRING) {
                //  if (na_len) { memcpy(ch, na_str, na_len); ch += na_len; }
                //} else if (quote) {
                //  QUOTE_FIELD;
                //} else {
                //  //tt0 = clock();
                //  memcpy(ch, CHAR(str), LENGTH(str));  // could have large fields. Doubt call overhead is much of an issue on small fields.
                //  ch += LENGTH(str);
                //  //tSTR += clock()-tt0;
                //}
                break;
              // default:
              // An uncovered type would have already thrown above when calculating maxLineLen earlier
              }
              *ch++ = col_sep;
              col_i++;
            }
          }
          ch--;  // backup onto the last col_sep after the last column
          memcpy(ch, row_sep, row_sep_len);  // replace it with the newline. TODO: replace memcpy call with eol1 eol2 --eolLen 
          ch += row_sep_len;
        }
      }
      #pragma omp ordered
      {
        if (f==-1) { 
          *ch='\0';  // standard C string end marker so Rprintf knows where to stop
          Rprintf(buffer);
          // despite being within '#pragma omp ordered' it seems Rprintf() even with a R_FlushConsole() too
          // still has some corruptions on Windows. Or it may be that test() uses a surrounding capture.output().
          // Anyway, fwrite.R now calls setDTthreads(1) when file="" (f==-1) which has solved 1729.3 and 1729.9.
        } else {
          WRITE(f, buffer, (int)(ch-buffer));
          // TODO: safe way to throw error from this thread if write fails (e.g. out disk space)
          //       { close(f); error("Error writing to file: %s", filename) };
          // Adding a progress bar here with Rprintf() or similar should be possible in theory, but see
          // the comment above about corruptions observed on windows. Could try and see. If the progress bar
          // corruputed then there's not too much harm. It could be delayed and flushed, or, best just have
          // one thread handle its update.
        }
        ch = buffer;
      }
    }
    free(buffer);
  }
  if (verbose) Rprintf("all %d threads done\n", nth);  // TO DO: report elapsed time since (clock()-t0)/NTH is only estimate
  if (f!=-1 && CLOSE(f)) error("Error closing file: %s", filename);
  return(R_NilValue);  // must always return SEXP from C level otherwise hang on Windows
}




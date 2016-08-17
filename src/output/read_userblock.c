#include <string.h>
#include <stdio.h>
extern char userblock_start;
extern char userblock_end;
extern char userblock_size;

long get_userblock_size_(void) 
{
   return (unsigned long)(&userblock_size);
}

long get_inifile_size(char* filename) 
{
   FILE* fp = fopen(filename, "rb");
   fseek(fp,0,SEEK_END);
   long length=ftell(fp);
   fclose(fp);
   return length;
}

void insert_userblock(char* filename, char* inifilename) 
{
   FILE* fp = fopen(filename, "rb+");
   rewind(fp);
   char* p = &userblock_start;
   while ( p != &userblock_end ) fputc(*p++, fp);
   fprintf(fp, "{[( INIFILE )]}\n");
   FILE* fini = fopen(inifilename, "rb");
   int c;
   do {
      c = fgetc (fini);
      if (c != EOF) fputc((char)c, fp);
   } while (c != EOF);
   fclose(fini);
   fprintf(fp, "{[( END USERBLOCK )]}\n");
   fclose(fp);
}


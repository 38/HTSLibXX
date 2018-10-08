#include <stdlib.h>
#include <stdio.h>
#include <api/BamReader.h>
#include <api/BamAlignment.h>


int main(int argc, char** argv)
{
	BamTools::BamReader br;
	br.Open(argv[1]);

	BamTools::BamAlignment al;

	while(br.GetNextAlignment(al))
	{
		uint32_t pos = al.Position;

		printf("%u\n", pos);
	}

	return 0;
}

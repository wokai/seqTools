/*
 * dna_fstream.h
 *
 *  Created on: 07.10.2013
 *      Author: wolfgang
 */

#ifndef DNA_FSTREAM_H_
#define DNA_FSTREAM_H_

# include <zlib.h>

/* + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
*
* 	DNA file stream:
*				Struct layer for reading text into char * buffer
*				from either uncompressed or compressed files.
*
* + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + */

/* + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
*
* 	DNA file stream:
*				Struct layer for reading text into char * buffer
*				from either uncompressed or compressed files.
*
* + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + */

static const int dfs_file_closed = 0;
static const int dfs_file_open   = 1;
static const int dfs_file_eof 	 = 2;

typedef struct dna_fstream
{
	unsigned type;
	int state;
	gzFile gz;

} dfStream;

static inline int dfs_isOpen(dfStream *dfs)	{ return dfs->state == dfs_file_open; }

static inline int dfs_stream_open(dfStream *dfs, const char *filename)
{
	dfs->gz = gzopen(filename, "rb");

		/*
		// Increase buffer size (from 8192 bytes)
		if(gzbuffer(dfs->gz,0xFFFFFF)<0) // 16.777.215
		{
			free(dfs->gz);
			dfs->gz=0;
		}
		*/


	if(dfs->gz)
		dfs->state = dfs_file_open;

	else
		dfs->state = dfs_file_closed;
	return dfs->state == dfs_file_closed;
}

static inline void dfs_stream_close(dfStream *dfs)
{
	if(dfs->state == dfs_file_open)
	{
		if(dfs->gz)
		{
			gzclose(dfs->gz);
			dfs->gz = 0;
			dfs->state = dfs_file_closed;
		}
		dfs->state = dfs_file_closed;
	}
}

static inline int dfs_stream_eof(dfStream *dfs)
{
	if(dfs->state == dfs_file_open)
	{
		if(dfs->gz)
			return gzeof(dfs->gz);
	}
	return dfs_file_eof;
}

static dfStream* dfs_stream_init(const char* filename)
{
	dfStream *dfs = calloc(sizeof(dfStream), 1);
	if(!dfs)
		return 0;
	dfs_stream_open(dfs, filename);
	return dfs;
}

static void dfs_destroy(dfStream *dfs)
{
	if(dfs)
	{
		if(dfs->state == dfs_file_open)
			dfs_stream_close(dfs);
		free(dfs);
	}
}

static inline size_t dfs_read(dfStream *dfs, char *dest, unsigned size)
{
	int res;
	size_t nchar;
	if(dfs_isOpen(dfs))
	{
		/*
		 * nchar<size: eof or error in both cases
		 */
		res=gzread(dfs->gz, dest, sizeof(char) * size);
		//printf("[dfs_read] Read size: %i\n",res);
		if(res < 0)
		{
			dfs_stream_close(dfs);
			return 0;
		}
		nchar= (size_t) res;

		if(dfs_stream_eof(dfs))
			dfs_stream_close(dfs);

		if( nchar < ((size_t) size))
			dfs_stream_close(dfs);

		return nchar;
	}
	return 0;
}
#endif /* DNA_FSTREAM_H_ */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include "include/graphtypes.h"
#include <assert.h>

#define MAX(x, y)   ((x) > (y) ? (x) : (y))

/* Entry index comparison function */
int poscmp(const void * p1, const void * p2)
{
    /* First compare contig name, then contig position */
    int chrresult = strncmp((*((pos_entry_t **)p1))->pos.contigname,
                            (*((pos_entry_t **)p2))->pos.contigname,
                            CONTIGNAME_LEN);
    if (chrresult == 0)
    {
        int64_t posdiff = (*((pos_entry_t **)p1))->pos.contigstart 
                - (*((pos_entry_t **)p2))->pos.contigstart;
        if (posdiff > 0) return 1;
        if (posdiff < 0) return -1;
        posdiff = (*((pos_entry_t **)p1))->pos.contigend 
                - (*((pos_entry_t **)p2))->pos.contigend;
        if (posdiff > 0) return 1;
        if (posdiff < 0) return -1;
        return 0;
    }
    
    return chrresult;
}

pos_entry_t * new_pos_entry(chr_pos_t * pos, uint64_t edgeindex)
{
    pos_entry_t * res = (pos_entry_t *) malloc(sizeof(pos_entry_t));
    if (pos != NULL)
    {
        res->pos = *pos;
    }
    res->edgeindex = edgeindex;
    
    return res;
}

map_entry_t * new_map_entry(pos_entry_t * entry)
{
    map_entry_t * res = (map_entry_t *) malloc(sizeof(map_entry_t));
    res->entry = entry;
    res->next = NULL;
    
    return res;
}

void realloc_edgemap(map_entry_t *** edgemap, uint64_t size,
        uint64_t * num_edges_allocated)
{
    /* Resize edgemap if necessary */
    if (size > *num_edges_allocated)
    {
        *edgemap = realloc(*edgemap, 
                *num_edges_allocated * 2 * sizeof(map_entry_t *));
        memset(*edgemap + *num_edges_allocated, 0,
                *num_edges_allocated * sizeof(map_entry_t *));
        *num_edges_allocated *= 2;
    }
}        

int has_node(map_entry_t * edgemap, map_entry_t * node)
{
    while (edgemap != NULL)
    {
        if (poscmp(&(edgemap->entry), &(node->entry)) == 0)
            return 1;
        edgemap = edgemap->next;
    }
    
    return 0;
}

void append_entry(map_entry_t ** edgemap, map_entry_t * node)
{
    if (!has_node(*edgemap, node))
    {
        node->next = *edgemap;
        *edgemap = node;
    }
}

/* Fills in a sorted repeat graph */
int main(int argc, char ** argv)
{
    /* Validate input, but not too much. */
    if (argc != 3)
    {
        fprintf(stderr, "usage: %s <graph_file> <contiglen>\n", argv[0]);
        return -1;
    }
    
    /* Open input */
    int input = open(argv[1], O_RDONLY);
    if (input < 0)
    {
        perror("opening input");
        return -1;
    }

    /* Get size of input file */
    struct stat inputstats;
    fstat(input, &(inputstats));
    off_t filesize = inputstats.st_size;
    
    /* Read input */
    uint8_t * filedata = (uint8_t *) malloc(filesize);
    if (filedata == NULL)
    {
        fprintf(stderr, "Out of memory\n");
        return -1;
    }
    off_t total = 0;
    while (total < filesize)
    {
        ssize_t count = read(input, filedata + total,
                filesize - total);
        if (count < 0)
        {
            perror("Reading input\n");
            return -1;
        }
        total += count;
    }
    
    /* Close input file */
    close(input);
    
    /* Calculate number of entries */
    uint64_t num_entries = filesize / sizeof(pos_entry_t);
    pos_entry_t * entries = (pos_entry_t *) filedata;

    /* Build map */
    uint64_t contiglen = atoll(argv[2]);
    uint64_t num_edges = 0;
    uint64_t num_edges_allocated = (num_entries / 2)+1;
    map_entry_t ** edgemap = (map_entry_t **) 
            malloc(num_edges_allocated * sizeof(map_entry_t *));
    assert(edgemap!=NULL);
    memset(edgemap, 0, num_edges_allocated * sizeof(map_entry_t *));
    for (uint64_t i = 0; i < num_entries; i++)
    {
        num_edges = MAX(num_edges, entries[i].edgeindex + 1);
        realloc_edgemap(&edgemap, num_edges, &num_edges_allocated);
        map_entry_t * me = new_map_entry(&(entries[i]));
        append_entry(&(edgemap[entries[i].edgeindex]), me);
    }

    /* Insert trivial "repeats" */
    char * lastcontig = "";
    uint64_t lastpos = 1;
    for (uint64_t i = 0; i < num_entries; i++)
    {
        pos_entry_t * cur_entry = &(entries[i]);
        
        /* If we've switched contigs */
        if (strcmp(lastcontig, cur_entry->pos.contigname))
        {
            /* If there's a gap at the start of this contig, */
            if (cur_entry->pos.contigstart != 1)
            {
                /* Create a new edge leading into this entry */
                chr_pos_t pos = { 0 };
                strcpy(pos.contigname, cur_entry->pos.contigname);
                pos.contigstart = 1;
                pos.contigend = cur_entry->pos.contigstart-1;
                pos.strand = '+';
                pos_entry_t * pe = new_pos_entry(&pos, num_edges++);
                realloc_edgemap(&edgemap, num_edges, &num_edges_allocated);
                map_entry_t * me = new_map_entry(pe);
                append_entry(&(edgemap[pe->edgeindex]), me);
            }
            /* cur_entry's memory isn't going anywhere... */
            lastcontig = cur_entry->pos.contigname;
            lastpos = cur_entry->pos.contigstart;
        }
        
        /* If there's a gap between this position and the last one, */
        if (cur_entry->pos.contigstart > lastpos)
        {
            /* Create a new edge to fill the gap */
            chr_pos_t pos = { 0 };
            strcpy(pos.contigname, cur_entry->pos.contigname);
            pos.contigstart = lastpos;
            pos.contigend = cur_entry->pos.contigstart-1;
            pos.strand = '+';
            pos_entry_t * pe = new_pos_entry(&pos, num_edges++);
            realloc_edgemap(&edgemap, num_edges, &num_edges_allocated);
            map_entry_t * me = new_map_entry(pe);
            append_entry(&(edgemap[pe->edgeindex]), me);

            lastpos = pos.contigend;
        }
        
        /* Check for funny business... */
        if (cur_entry->pos.contigstart < lastpos)
        {
            fprintf(stderr, "Funny edge at %s:%lu->%lu\n",
                    cur_entry->pos.contigname,
                    cur_entry->pos.contigstart,
                    cur_entry->pos.contigend);
        }
        
        /* Advance lastpos */
        lastpos = cur_entry->pos.contigend+1;
    }
    
    /* Fill in up to contig len */
    if (lastpos < contiglen)
    {
        chr_pos_t pos = { 0 };
        strcpy(pos.contigname, lastcontig);
        pos.contigstart = lastpos;
        pos.contigend = contiglen;
        pos.strand = '+';
        pos_entry_t * pe = new_pos_entry(&pos, num_edges++);
        realloc_edgemap(&edgemap, num_edges, &num_edges_allocated);
        map_entry_t * me = new_map_entry(pe);
        append_entry(&(edgemap[pe->edgeindex]), me);
    }
    
    /* Sort the graph */
    uint64_t index_size = 0;
    for (uint64_t i = 0; i < num_edges; i++)
    {
        map_entry_t * me = edgemap[i];
        while (me != NULL)
        {
            index_size++;
            me = me->next;
        }
    }
    
    pos_entry_t ** pos_index = (pos_entry_t **) 
            malloc(index_size * sizeof(pos_entry_t *));
    index_size = 0;
    for (uint64_t i = 0; i < num_edges; i++)
    {
        map_entry_t * me = edgemap[i];
        while (me != NULL)
        {
            pos_index[index_size++] = me->entry;
            me = me->next;
        }
    }
    qsort(pos_index, index_size, sizeof(pos_entry_t *), poscmp);
  
    /* Compute prefix sum on edges */
    uint64_t * ps_edgemap = (uint64_t *) malloc(num_edges * sizeof(uint64_t));
    uint64_t cur_edge = 0;
    for (uint64_t i = 0; i < num_edges; i++)
    {
        /* Exclusive */
        ps_edgemap[i] = cur_edge;
        if (edgemap[i] != NULL)
            cur_edge++;
    }
    
    /* Write out graph */
    int output = open(argv[1], O_WRONLY | O_TRUNC);
    if (output < 0)
    {
        perror("opening output");
        return -1;
    }
    
    for (uint64_t i = 0; i < index_size; i++)
    {
        /* Adjust edge index */
        pos_index[i]->edgeindex = ps_edgemap[pos_index[i]->edgeindex];
        
        /* Write out to output */
        write(output, pos_index[i], sizeof(pos_entry_t));
    }
    
    /* Close output */
    close(output);
            
    return 0;
}

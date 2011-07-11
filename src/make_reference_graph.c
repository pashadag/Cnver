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
#include <assert.h>
#include "include/graphtypes.h"

#define MIN_EDGE_LEN    2

uint64_t INVALID_EDGE=-1;

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

/* Intersection function */
int posintersect(const chr_pos_t * p1, const chr_pos_t * p2)
{
    if (strncmp(p1->contigname, p2->contigname, CONTIGNAME_LEN) != 0)
        return 0;
    return (p1->contigend >= p2->contigstart);
}

/* Utility functions */
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

#ifdef DUMP
FILE * f_out=NULL;
void print_pos(chr_pos_t *p) {
	fprintf(f_out,"POS: %llu %llu %c\n",p->contigstart,p->contigend,p->strand);
}


void print_edge(map_entry_t ** edgemap, uint64_t edge) {
	if (edge==INVALID_EDGE) {
		fprintf(f_out,"INVALID EDGE\n");
		return;
	}
	
	fprintf(f_out,"edge: %llu\n",edge);
	map_entry_t * me =edgemap[edge];
	while (me!=NULL) {
		fprintf(f_out,"\t%llu\t%llu\t%c\n",me->entry->pos.contigstart, me->entry->pos.contigend, me->entry->pos.strand);
		me=me->next;
	}

}
void write_out(uint64_t num_edges, map_entry_t ** edgemap,uint64_t index_size,pos_entry_t ** pos_index, uint64_t iteration) {
    /* Write out graph */
    char filename[200];
    sprintf(filename,"mkgraph.dump.%0.2d",iteration);
    if (f_out==NULL) {
    	f_out=fopen(filename,"w");
    }
    if (f_out==NULL) {
	printf("ERROR IN LOGGING\n");
	exit(1);
    }

    
    fprintf(f_out,"Pos_index:\n");
    for (uint64_t i =0; i< index_size; i++) {
 	print_pos(&(pos_index[i]->pos));
    }

    fprintf(f_out,"Edge_map:\n");
    for (uint64_t i = 0; i < num_edges; i++)
    {
	map_entry_t * t=edgemap[i];
	while (t!=NULL) {	
	   fprintf(f_out,"E: %8llu\tS:%llu\tE:%llu\t%c\n",t->entry->edgeindex,t->entry->pos.contigstart,t->entry->pos.contigend,t->entry->pos.strand);       
	   t=t->next;
	}
    }
}


#endif
uint64_t split_edge(map_entry_t *** pointer_to_edgemap,uint64_t old_edge,uint64_t * num_edges, uint64_t * num_edges_allocated,uint64_t offset,char strand) {
		/*
			old_edge
		|-----------------------------|
		|---------------  --------------|
		|---------------|         |--------------|
		    left_edge (old_edge)   right_edge (new edge)
		*/
	    map_entry_t ** edgemap = (*pointer_to_edgemap);
            if (edgemap[old_edge]==NULL) {
		printf("EDGE DOES NOT EXIST!");
		exit(1);
	    }

	    chr_pos_t any_pos = edgemap[old_edge]->entry->pos;
	    uint64_t old_edge_len = any_pos.contigend - any_pos.contigstart+1;
	    uint64_t right_edge_len = old_edge_len - offset;
	    uint64_t right_edge = (*num_edges)++;
	    realloc_edgemap(pointer_to_edgemap, *num_edges, num_edges_allocated);
	    //assert((*num_edges_allocated)>(*num_edges));
	    edgemap = (*pointer_to_edgemap);
	    //assert((*num_edges_allocated)>right_edge);
	    map_entry_t * me = edgemap[old_edge];

	    while (me != NULL)
	    {
		    //The current pos entry
		    chr_pos_t * current_pos = &(me->entry->pos);
		    pos_entry_t * new_pos = new_pos_entry(current_pos,right_edge);
#ifdef DUMP
			fprintf(f_out,"Splitting pos\n");
			print_pos(current_pos);
#endif
		    if (current_pos->strand == strand) {
			new_pos->pos.contigstart = current_pos->contigstart + offset;
		    } else {
			new_pos->pos.contigend = current_pos->contigend - offset;
		    }
		    map_entry_t * copy = new_map_entry(new_pos);
		    append_entry(&(edgemap[right_edge]), copy);
		
		    if (current_pos->strand == strand)
		    {
		       current_pos->contigend = current_pos->contigstart + offset - 1;
		    } else {
		       current_pos->contigstart = current_pos->contigend - offset + 1;
		    }
#ifdef DUMP
			fprintf(f_out,"TO:\n");
			print_pos(current_pos);
			print_pos(&(new_pos->pos));
#endif
		    /* Move on to next entry */
		    me = me->next;
	    }

}

void merge_edges(map_entry_t ** edgemap,uint64_t from_edge,uint64_t to_edge,uint64_t num_edges,int flip) {
    if (from_edge == INVALID_EDGE || to_edge == INVALID_EDGE ) {
	return;
    }
    if (from_edge==to_edge) {

	printf("Merging the same edge???\n");
	exit(1);
    }
    map_entry_t * mover = edgemap[from_edge];
    while (mover != NULL)
    {
	/* Don't move already existing nodes */
	if (has_node(edgemap[to_edge], mover))
	{
	    /* Mark node as unused */
	    mover->entry->edgeindex = INVALID_EDGE;
	    map_entry_t * tmp = mover->next;
	    free(mover);
	    mover = tmp;
	}
	else
	{
	    if (flip==1) {
		if (mover->entry->pos.strand=='+') {
			mover->entry->pos.strand='-';
		} else {
			mover->entry->pos.strand='+';
		}	
	    }
	    mover->entry->edgeindex = to_edge;
	    map_entry_t * tmp = mover->next;
	    mover->next = edgemap[to_edge];
	    edgemap[to_edge] = mover;
	    mover = tmp;
	}
    }

    edgemap[from_edge]=NULL;

}


uint64_t reindex_nodes(map_entry_t ** edgemap, pos_entry_t *** pos_index, uint64_t num_edges) {
        /* Re-index nodes */
        uint64_t index_size = 0;
        for (uint64_t i = 0; i < num_edges; i++)
        {
            map_entry_t * me = edgemap[i];
            while (me != NULL)
            {
		assert(me->entry->edgeindex!=INVALID_EDGE);
		//printf("ADDING\n");
		//print_edge(edgemap,i);
                index_size++;
                me = me->next;
            }
        }
	//printf("INDEX SIZE IS %llu\n",index_size);
        *pos_index = realloc(*pos_index, index_size * sizeof(pos_entry_t **));
        index_size = 0;
        for (uint64_t i = 0; i < num_edges; i++)
        {
            map_entry_t * me = edgemap[i];
            while (me != NULL)
            {
                (*pos_index)[index_size++] = me->entry;
                me = me->next;
            }
        }
        qsort(*pos_index, index_size, sizeof(pos_entry_t *), poscmp);
	//printf("INDEX SIZE IS %llu\n",index_size);
	return index_size;
}


int merge_duplicates(pos_entry_t ** pos_index,map_entry_t ** edgemap,uint64_t index_size, uint64_t num_edges) { 
    int removed=0;
    /* Remove duplicate edges */
    for (int i = 1; i < index_size; i++)
    {
        pos_entry_t * prev_node = pos_index[i-1];
        pos_entry_t * cur_node = pos_index[i];
        
        /* Skip unused entries */
        if (cur_node->edgeindex == num_edges)
            continue;
            
        uint64_t from_edge = prev_node->edgeindex;
        uint64_t to_edge = cur_node->edgeindex;
        /* Do pairwise comparison */
        if (poscmp(&prev_node, &cur_node) == 0 && (from_edge < INVALID_EDGE && to_edge < INVALID_EDGE) )
        {
		
#ifdef DUMP
    fprintf(f_out,"Merging edge %d with %d:\n",to_edge,from_edge);
    fprintf(f_out,"to_edge\n");
    print_pos(&(cur_node->pos));
    print_edge(edgemap,to_edge);
    fprintf(f_out,"from_edge\n");
    print_edge(edgemap,from_edge);
    print_pos(&(prev_node->pos));
    fprintf(f_out,"done merge\n");
#endif
	    //printf("FOUND SAME:\n");
	    //print_edge(edgemap,prev_node->edgeindex);
	    //print_edge(edgemap,cur_node->edgeindex);
            /* Move current node's edgemates to previous node's edge */
	    assert(to_edge<num_edges);
	    assert(from_edge<num_edges);
	    /* Check if we need to flip the new coordinates */
	    int flip=0;
	    if (prev_node->pos.strand!=cur_node->pos.strand) {
		flip=1;
	    }
	    removed++;
	    merge_edges(edgemap,from_edge,to_edge,num_edges,flip);
        }
   }
   return removed;
}


/* Takes, as input, repeat graph entries. Produces, as output,
    a sorted repeat graph */
int main(int argc, char ** argv)
{
	//printf("INVALID EDGE IS %llu\n",INVALID_EDGE);
    /* Validate input, but not too much. */
    if (argc != 3)
    {
        fprintf(stderr, "usage: %s <bin_file> <output_file>\n", argv[0]);
        return -1;
    }
#ifdef DUMP

	printf("DUMP !! \n");
#endif
    
    /* Open input and output files */
    int input = open(argv[1], O_RDONLY);
    int output = open(argv[2], O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
    if ((input < 0) || (output < 0))
    {
        perror("opening input/output");
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

    /* Build edge map */
    /* Number of edges at this point is just half the number of entries */
    uint64_t num_edges = num_entries / 2;
    map_entry_t ** edgemap = 
            (map_entry_t **) malloc(num_edges * sizeof(map_entry_t *));
    for (int i = 0; i < num_edges; i++)
    {
        edgemap[i] = new_map_entry(entries + (i * 2));
        edgemap[i]->next = new_map_entry(entries + (i * 2) + 1);
    }
    
    /* Make position index */
    pos_entry_t ** pos_index = 
            (pos_entry_t **) malloc(num_entries * sizeof(pos_entry_t *));
    for (int i = 0; i < num_entries; i++)
        pos_index[i] = entries + i;
    qsort(pos_index, num_entries, sizeof(pos_entry_t *), poscmp);



    /* Remove duplicate edges */
    merge_duplicates(pos_index,edgemap,num_entries,num_edges);

    
    uint64_t num_edges_allocated = num_edges;
    int convergence;
    int iterations = 0;
    uint64_t index_size = 0;
    do {
	//printf("ITERATION %d\n",iterations);	
        convergence = 1;
        iterations++;
#if 0        
        /* Cull short edges */
        for (uint64_t i = 0; i < num_edges; i++)
        {
            map_entry_t * me = edgemap[i];
            if ((me != NULL) &&
                    ((me->entry->pos.contigend - 
                        me->entry->pos.contigstart + 1) < MIN_EDGE_LEN))
            {
                edgemap[i] = NULL;
            }
        }
#endif   


	/*

		This is a 4 step process
		1) Go over all edges and allocate space only for used positions, fill this space with pointers
		to each of the positions, sort the pointers in place
		2) If found two duplicate positions, merge the two edges they correspond to together
		3) Repeat step 1, since some positions may have been removed (during merging)
		4) Split overlapping edges
			a) If have this condition
			
			|---------------------------- ...  (Prev)
			      |---------------------- ...  (Cur)

			Then split and do this

			|----||---------------------- ... (Prev is first one) (New edge positions are not in this index
					but will be on next iteration)


			b)  If have this condition

			... --------------|		(Prev)
			... ---------------------|	(Cur)

			Then split and do this

			... --------------|		(Prev)
			... --------------||-----|	(Cur is first one) (New edge positions are no in this index but will
					be on next iteration)

			NOTE:
				Any overlapping edges may pass through a and/or b but then will pass through 2 and get merged
				This should account for all cases


	*/


	index_size=reindex_nodes(edgemap,&pos_index,num_edges); 
#ifdef DUMP
	if (f_out!=NULL){ 
		fclose(f_out);
		f_out=NULL;
	}
	write_out(num_edges,edgemap,index_size,pos_index,iterations);
#endif 
        /* Remove duplicate edges */
        int removed=merge_duplicates(pos_index,edgemap,index_size,num_edges);
	index_size=reindex_nodes(edgemap,&pos_index,num_edges); 

        //fprintf(stderr, "[%lu]: %lu nodes, %lu edges, r %d\n", iterations, index_size, num_edges,removed);

        /* Break apart coincident edges */
        for (uint64_t i = 1; i < index_size; i++)
        {
            pos_entry_t * prev_entry = pos_index[i-1];
            pos_entry_t * cur_entry = pos_index[i];

            chr_pos_t prev_pos = prev_entry->pos;
            chr_pos_t cur_pos = cur_entry->pos;

	    if  ((poscmp(&prev_entry, &cur_entry) < 0)  && (posintersect(&prev_pos, &cur_pos))) {
                /* Get differences in edge positions */
                int64_t cs_diff = cur_pos.contigstart - prev_pos.contigstart;
                int64_t ce_diff = cur_pos.contigend - prev_pos.contigend;
                uint64_t prev_len = prev_pos.contigend - prev_pos.contigstart + 1;
                uint64_t cur_len = cur_pos.contigend - cur_pos.contigstart + 1;


		uint64_t prev_index = prev_entry->edgeindex;
		uint64_t cur_index = cur_entry->edgeindex;
		assert(prev_index<num_edges);
		assert(cur_index<num_edges);


#ifdef DUMP
		    fprintf(f_out,"\nConsidering,\n");
		    print_pos(&cur_pos);
		    print_pos(&prev_pos);
#endif

                /* Different types of intersection: */
                /* Complete overlap, coincident start */
                if (cs_diff != 0)
                {
		    /*  
			current_start > previous_start
		    */
		    split_edge(&edgemap,prev_entry->edgeindex,&num_edges,&num_edges_allocated,cs_diff,prev_pos.strand);
                } else if (ce_diff > 0) {
		    /*
			current_end > previous_end 
		    */
		    split_edge(&edgemap,cur_entry->edgeindex,&num_edges,&num_edges_allocated,cur_len-ce_diff,cur_pos.strand);
                }
                
                convergence = 0;
            }

        }
    } while (!convergence);
   
#ifdef DUMP
	if (f_out!=NULL) {
		fclose(f_out);
		f_out=NULL;
	}

#endif

 
    fprintf(stderr, "Converged after %lu iterations\n", iterations);
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
    
    fprintf(stderr, "num_edges: %lu\nnum_nodes: %lu\n",
            cur_edge, index_size);


    /* Write out graph */
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

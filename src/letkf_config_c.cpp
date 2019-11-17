// Copyright 2019-2019 Travis Sluka
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


/**
 * \cond internal
 * Interfacing between the fortran and C++ yaml library
 */

#include "yaml.h"
#include <iostream>

extern "C" {

void letkf_yaml_get_value(yaml_node_t *&node, char * str, int &str_len){
        // make sure this is a scalar node
        if (node->type != YAML_SCALAR_NODE) {
                str_len = -1;
                return;
        }

        //copy the string
        strcpy(str,(char*)node->data.scalar.value);
        str_len = node->data.scalar.length;
}

int letkf_yaml_init( yaml_document_t *&yaml_doc, yaml_node_t *&root_node,
                     char *filename ){
        yaml_parser_t parser;
        FILE *f;
        int i;

        //std::cout <<"DBG: loading file: "<<filename<<std::endl;
        // read in and parse the yaml document
        f=NULL;
        yaml_doc = new yaml_document_t;
        f=fopen(filename, "r");
        if (f == NULL) return 1;
        if (yaml_parser_initialize(&parser) == 0) return 1;
        yaml_parser_set_input_file(&parser, f);
        if (yaml_parser_load(&parser, yaml_doc) == 0) return 2;
        yaml_parser_delete(&parser);
        root_node = yaml_document_get_root_node(yaml_doc);
        fclose(f);

        // read in other '#include' files
        int stack_len = 0;

        yaml_node_t *stack[1024], *ptr;
        stack[stack_len] = root_node;
        stack_len++;
        while (stack_len > 0) {
                ptr = stack[stack_len-1];
                stack_len--;

                // if a scalar, skip
                if (ptr->type == YAML_SCALAR_NODE) continue;

                // if a sequence, add each child to the stack
                if (ptr->type == YAML_SEQUENCE_NODE) {
                        yaml_node_item_t *ptr2 = ptr->data.sequence.items.start;
                        while(ptr2 != ptr->data.sequence.items.top) {
                                stack[stack_len] = yaml_document_get_node(yaml_doc, *ptr2);
                                stack_len++;
                                ptr2++;
                        }
                }

                // if a mapping, add each child to stack, OR
                if (ptr->type == YAML_MAPPING_NODE) {
                        yaml_node_pair_t *node_kv;
                        yaml_node_t *node_k, *node_v;
                        node_kv = ptr->data.mapping.pairs.start;
                        while (node_kv != ptr->data.mapping.pairs.top) {
                                node_k = yaml_document_get_node(yaml_doc, node_kv->key);
                                node_v = yaml_document_get_node(yaml_doc, node_kv->value);
                                if (node_k->type != YAML_SCALAR_NODE) {
                                        stack[stack_len]=node_v;
                                        stack_len++;
                                }
                                else {
                                        // all keys should be scalar nodes... but just to be safe
                                        //std::cout <<"searching "<<node_k->data.scalar.value<<std::endl;
                                        if (strcmp((char*)node_k->data.scalar.value, "#import") == 0) {
                                          std::cout <<"ERROR: importing files not yet implemented"<<std::endl;
                                          return 1;
                                                char str[1024];
                                                int str_len;
                                                yaml_document_t *doc2 = new yaml_document_t;
                                                yaml_node_t *root2;
                                                letkf_yaml_get_value(node_v, str, str_len);
                                                letkf_yaml_init(doc2, root2, str);
                                        }
                                }
                                node_kv++;
                        }
                }
        }
        return 0;
}



int letkf_yaml_get_len(yaml_document_t *&doc, yaml_node_t *&node){
        int count = 0;
        if (node->type == YAML_SEQUENCE_NODE) {
                yaml_node_item_t *idx=node->data.sequence.items.start;
                while(idx != node->data.sequence.items.top) {
                        count++;
                        idx++;
                }
        }
        else if(node->type == YAML_MAPPING_NODE) {
                yaml_node_pair_t *idx=node->data.mapping.pairs.start;
                while(idx != node->data.mapping.pairs.top) {
                        count++;
                        idx++;
                }
        }
        return count;
}



yaml_node_t * letkf_yaml_get_child_name(yaml_document_t *&doc,
                                        yaml_node_t *&node, char *key){
        yaml_node_pair_t *node_kv;
        yaml_node_t *node_k, *node_v;

        // abort if this is not a mapping
        if(node->type == YAML_SCALAR_NODE) return NULL;
        if(node->type == YAML_SEQUENCE_NODE) return NULL;

        //find the matching key
        node_kv = node->data.mapping.pairs.start;
        while (node_kv != node->data.mapping.pairs.top) {
                node_k = yaml_document_get_node(doc,  node_kv->key);
                node_v = yaml_document_get_node(doc,  node_kv->value);
                if (node_k == NULL) return NULL;
                if (node_v == NULL) return NULL;
                if (node_k->type == YAML_SCALAR_NODE) {
                        // all keys should be scalar nodes... but just to be safe
                        if (strcmp((char*)node_k->data.scalar.value, key) == 0) {
                                return node_v;
                        }
                }
                node_kv++;
        }
        return NULL;
}



yaml_node_t * letkf_yaml_get_child_idx(yaml_document_t *&doc,
                                       yaml_node_t *&node, int &idx,
                                       char * name_str, int &name_len){

        int maxlen = letkf_yaml_get_len(doc, node);
        if (idx > maxlen) return NULL;

        yaml_node_t *res_node = NULL;
        if (node->type == YAML_SEQUENCE_NODE) {
                yaml_node_item_t *n2 = node->data.sequence.items.start;
                for(int i=1; i != idx; i++) n2++;
                res_node = yaml_document_get_node(doc, *n2);
                name_len = -1;
        }
        else if(node->type == YAML_MAPPING_NODE) {
                yaml_node_pair_t *n2 = node->data.mapping.pairs.start;
                for(int i=1; i != idx; i++) n2++;
                res_node = yaml_document_get_node(doc, n2->value);
                yaml_node_t * node_k = yaml_document_get_node(doc, n2->key);
                //TODO make sure this is a scalar
                name_len = node_k->data.scalar.length;
                strcpy(name_str,(char*)node_k->data.scalar.value);
        }
        return res_node;

}

}

/**
 * \endcond
 **/

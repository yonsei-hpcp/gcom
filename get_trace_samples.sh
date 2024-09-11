#!/bin/bash
if [ ! -d trace_samples ]; then
	wget https://hpcp.yonsei.ac.kr/~jounghoo/GCoM/trace_samples.tar.gz
	tar xzvf trace_samples.tar.gz
	rm trace_samples.tar.gz
fi

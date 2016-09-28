#! /usr/bin/bash

hostname
echo "hogehoge" > /mnt/scratch/test.$$
hostname >> /mnt/scratch/test.$$
echo "fugafuga" >>  /mnt/scratch/test.$$
sleep 10

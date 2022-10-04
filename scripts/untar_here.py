#!/usr/bin/env python

"""
untar_here.py
-------------
Bundled script to un-tar the eniric data downloaded.

Uses the tarfile module to extract the data.

"""
import argparse
import sys
import tarfile


def _parser():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(description="Extract from a tar file.")

    parser.add_argument("filename", help="File to untar.", type=str, default="")
    return parser.parse_args()


if __name__ == "__main__":
    filename = _parser().filename

    with tarfile.open(filename, "r") as tar:
        
        import os
        
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar)
    sys.exit(0)

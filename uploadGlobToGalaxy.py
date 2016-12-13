#!/usr/bin/env python
from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.objects import GalaxyInstance as GI
import argparse
import sys
import os
import glob


def getLibID(gi, libName):
    """
    Given a library name, like "foo bar", return the ID for the first library matching that
    """
    lib = gi.libraries.get_libraries(name=libName)
    if not lib or len(lib) == 0:
        libs = gi.libraries.get_libraries()
        # Handle differences in capitalization
        for l in libs:
            if l["name"].lower() == libName.lower():
                return l["id"]
        raise RuntimeError("No library named {}".format(libName))
    return lib[0]["id"]


def getFolderID(gi, libID, path):
    """
    Given a library ID (lib) and a path (e.g., "/foo/bar/sniggly3"), return the folder ID for sniggly3.
    
    If the path doesn't exist (in part or in total), then create it.
    """
    # Does the path already exist?
    folders = gi.libraries.get_folders(libID, name=path)
    if folders is not None and len(folders) > 0:
        return folders[0]["id"]
    
    # Get the closest base folder
    longest = gi.libraries.get_folders(libID, name="/")[0]
    folders = gi.libraries.get_folders(libID)
    for folder in folders:
        if path.startswith(folder["name"]):
            # Look for the longest pathname overlap. The next character in path MUST be "/",
            # since otherwise adding "/a/b/c2" when "/a/b/c" exists would result in "/a/b/c/2"!
            if len(folder["name"]) > len(longest["name"]) and path[len(folder["name"])] == "/":
                longest = folder
    
    # shorten the path name if relevant
    idx = len(longest["name"])
    pathLeft = path[idx:]
    if pathLeft.startswith("/"):
        pathLeft = pathLeft[1:]

    for fName in pathLeft.split("/"):
        gi.libraries.create_folder(libID, fName, base_folder_id=longest["id"]) # returns None
        if longest["name"] != "/":
            newFName = "{}/{}".format(longest["name"], fName)
        else:
            newFName = "/{}".format(fName)
        longest = gi.libraries.get_folders(libID, name=newFName)[0]
    
    return longest["id"]


def addFileToLibraryFolder(gi, libID, folderID, fileName, file_type='auto', dbkey='?', link=True, roles=''):
    """
    Link/copy "fname" into the library and folder specified by libID and folderID. These MUST exist.
    
    file_type, dbkey, and roles are pass through to upload_from_galaxy_filesystem().
    
    link must be True (default) or False. If it's False then files are copied in.
    
    This returns a dictionary with keys: name, url, and id (or presumably None on error).
    """
    if link == True:
        link_data_only = 'link_to_files'
    else:
        link_data_only = 'copy_files'
    
    rv = gi.libraries.upload_from_galaxy_filesystem(libID,
                                                    fileName,
                                                    folder_id=folderID,
                                                    file_type=file_type,
                                                    dbkey=dbkey,
                                                    link_data_only=link_data_only,
                                                    roles=roles)
    return rv


def addFileToLibrary(gi, libraryName, path, fileName, file_type='auto', dbkey='?', link=True, roles=''):
    """
    Add fileName to path in libraryName. gi is a GalaxyInstance
    
    The other parameters are passed to upload_from_galaxy_filesystem()
    """
    libID = getLibID(gi, libraryName)
    folderID = getFolderID(gi, libID, path)
    dataset = addFileToLibraryFolder(gi, libID, folderID, fileName, file_type=file_type, dbkey=dbkey, link=link, roles='')
    if not dataset:
        raise RuntimeError("Error adding '{}' to '{}'".format(fileName, path))
    return dataset


def getFileType(fName):
    """
    If the file name ends with .fastq.gz then return 'fastqsanger'. Otherwise, return 'auto'.
    """
    if fName.endswith(".fastq.gz") or fName.endswith(".fq.gz"):
        return "fastqsanger"
    return "auto"


def checkExists(datasets, folderID, fName):
    """
    Return true if there's a file named fName in a folder with folderID in side the library with ID libID. Otherwise, return False
    For fastq.gz files, the file name in Galaxy might be lacking the .gz extension when we check...
    """
    for dataset in datasets:
        if dataset.wrapped["folder_id"] == folderID:
            if dataset.name == os.path.basename(fName):
                return True
            if fName.endswith(".fastq.gz") and os.path.basename(fName)[:-3] == dataset.name:
                return True
    return False


def getFiles(d):
    """
    Given a directory, return a list of all files in all subdirectories
    """
    files = []
    for (dpath, dnames, filenames) in os.walk(d):
        for fname in filenames:
            files.append("{}/{}".format(dpath, fname))
    return files


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", "-f", action="store_true", default=False, help="If specified, files already present will be added again.")
    parser.add_argument("--basePath", help="This is removed from each file name and the remainder is appended to 'libName' before adding each file. This is convenient with globbing multiple folders where the folder structure should be maintained.")
    parser.add_argument("libName", help="The library name (under shared-data -> data libraries")
    parser.add_argument("path", default="/", help="The path of folders to put the file/files into.")
    parser.add_argument("fglobs", help="file name(s), which can also be globs (e.g., '*.fastq.gz')", nargs="+")
    return parser


def main(args=None):
    args = parseArgs(args).parse_args()

    # Read in the userKey and set the server
    f = open("{}/.galaxy_key".format(os.path.expanduser("~")))
    userKey = f.readline().strip()
    f.close()
    url = "http://galaxy.ie-freiburg.mpg.de"

    gi = GalaxyInstance(url=url, key=userKey)
    # There's no way to get datasets in a library without using the wrapper clients
    gi2 = GI(url=url, api_key=userKey)

    # memoize the library ID
    libID = getLibID(gi, args.libName)

    # Use the correct library name
    _ = gi.libraries.get_libraries(library_id=libID)[0]["name"]
    l = gi2.libraries.list(name=_)[0]
    currentDatasets = l.get_datasets()

    # path needs to start with a "/"
    if not args.path.startswith("/"):
        args.path = "/{}".format(args.path)

    # memoize the folder ID
    folderID = getFolderID(gi, libID, args.path)

    for foo in args.fglobs:
        # handle any globs
        globs = glob.glob(foo)

        for g in globs:
            # If we have a directory, make a list of files
            # If we have a file, put it in a list
            if os.path.isdir(g):
                fileList = getFiles(g)
                if not args.basePath:
                    args.basePath = g
            else:
                fileList = [g]

            # Trim "/" from args.basePath
            if args.basePath:
                args.basePath = args.basePath.rstrip("/")

            for fName in fileList:
                # iterate over the files
                if args.basePath:
                    # If we have a base path, then deal with it.
                    bname = os.path.dirname(fName)
                    assert(fName.startswith(args.basePath))
                    newPath = args.path + bname[len(args.basePath):]
                    # if args.path == "/", though, we end up with "//something..."
                    if newPath.startswith("//"):
                        newPath = newPath[1:]
                    folderID = getFolderID(gi, libID, newPath)
    
                # If the file already exists then don't add multiple copies unless --force is used.
                if not args.force:
                    if checkExists(currentDatasets, folderID, fName):
                        sys.stderr.write("Skipping {}, already added\n".format(os.path.basename(fName)))
                        continue
                f = addFileToLibraryFolder(gi, libID, folderID, fName, file_type=getFileType(fName))
                if not f:
                    sys.exit("There was an error while adding {}".format(fName))


if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(args)

#!/usr/bin/env python3
import argparse
import sys
import smtplib
import os
import smtplib
from email.mime.text import MIMEText


def loadUserDictionary():
    d = dict()
    f = open("/home/pipegrp/parkourUsers.txt")
    for line in f:
        cols = line.split("\t")
        d[cols[1]] = [cols[0], cols[2]]
    f.close()
    return d


def fetchEmailAndName(project, userDictionary):
    lastName = project.split("_")[-2]
    if lastName not in userDictionary:
        return None, None
    else:
        return userDictionary[lastName]

def getProjectIDs(projects):
    IDs = []
    for p in projects:
        # Sanity check
        assert(p.startswith("Project_"))
        IDs.append(p.split("_")[1])
    if len(IDs) == 1:
        return IDs[0]
    return " and ".join([", ".join(IDs[:-1]), IDs[-1]])


def getFlowCell():
    return os.path.split(os.getcwd())[-1]


def main(args):
    parser = argparse.ArgumentParser(description="Send an email to one or more users about a project or project being finished. This must be run in the output directory of the demultiplexing pipeline.")
    parser.add_argument("--notGood", action="store_true", help="If specified, do NOT say that the sequencing quality was good.")
    parser.add_argument("--analysis", action="store_true", help="If specified, the BigRedButton did something with these projects.")
    parser.add_argument("--cc", nargs="+", help="One or more addresses to CC.")
    parser.add_argument("--comment", help="Either comment that will be included as its own paragraph (ensure you quote the whole thing!) or the path to a file containing such a comment.")
    parser.add_argument("--noGalaxy", action="store_true", help="Do NOT say that files are in Galaxy.")
    parser.add_argument("--fromPerson", help="The name of the person sending the email. Note that the contents of ~/.signature are also used! (Default: Devon)", default="Devon")
    parser.add_argument("--fromEmail", help="The email address of the person sending this. Note that they receive a copy as BCC!", default="ryan@ie-freiburg.mpg.de")
    parser.add_argument("project", nargs="+", help="One or more project directories. Only the user on the first will receive an email!")
    args = parser.parse_args(args)

    userDictionary = loadUserDictionary()
    firstName, email = fetchEmailAndName(args.project[0], userDictionary)
    if not firstName:
        sys.exit("That user isn't known!\n")

    content = """Hi {},

Your sequencing samples for project""".format(firstName)

    if len(args.project) > 1:
        content += "s"
    content += " {} are finished and the results are now available in your group's sequencing data directory".format(getProjectIDs(args.project))

    if not args.noGalaxy:
        content += " and Galaxy"
    content += " under the {} folder.".format(getFlowCell())

    if not args.notGood:
        content += " The overall sequencing quality for these samples was good."
    if args.analysis:
        content += " An automated partial analysis (https://www.biorxiv.org/content/early/2018/09/04/407312) of these samples is available in the same location. If you would like completed analysis of these samples in a semi-automated fashion, please request that using our online portal: http://snakequest.ie-freiburg.mpg.de .\n"

    if args.comment:
        content += "\n"
        if os.path.exists(args.comment):
            content += open(args.comment).read()
        else:
            content += args.comment
        content += "\n"

    content += "\n\nPlease let me know if you have any questions,\n{}\n".format(args.fromPerson)

    # Add a .signature
    if os.path.exists("/home/pipegrp/.signature"):
        content += "\n--\n"
        content += open("/home/pipegrp/.signature").read()

    # Send the Email
    msg = MIMEText(content)
    msg['Subject'] = "Sequencing samples ready"
    msg['From'] = args.fromEmail
    msg['To'] = email
    if args.cc:
        msg['Cc'] = ", ".join(args.cc)
    msg['Bcc'] = args.fromEmail

    s = smtplib.SMTP("mail.ie-freiburg.mpg.de")
    s.send_message(msg)
    s.quit()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("--help")
    main(sys.argv[1:])

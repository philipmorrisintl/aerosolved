#!/bin/bash

# Generates aerosolModel/aerosolModelGitInfo.H with the current git state so it
# can be printed to the log by aerosolModel.C.
#
# Portable-across-GNU-BSD:
#   * staleness is tracked with a commit sidecar (.aerosolvedGitInfoCommit)
#     rather than `find -mmin` / `touch -d`, which are GNU-only
#   * body lines are emitted with `printf '%s\n'` instead of `echo "...."`;
#     this avoids the xpg_echo / POSIX `echo` backslash-interpretation
#     ambiguity that produced over-escaped \" sequences on some distros
#     (see IMPROVEMENTS.txt TIER-1 #A-7 for the diagnosis)
#
# The .H is regenerated ONLY when the current git state (branch + commit +
# working-tree diff stat signature) differs from the last recorded signature.
# When git is absent or not a repository, it emits a single no-git stub line.

AEROSOLVEDROOT=../../
GITINFOFILE=aerosolModel/aerosolModelGitInfo.H
COMMITFILE=aerosolModel/.aerosolvedGitInfoCommit
TEMPGITINFOFILE=tempGitInfoFile.txt

GITBRANCH=""
GITCOMMIT=""

# Decide whether the git info must be regenerated.
# We key on a signature that contains branch + commit + diffstat; this changes
# only when something the .H actually prints would change, so it stays
# idempotent across repeat builds.
NEEDS=1

if [ -x "$(command -v git)" ] && [ -d "$AEROSOLVEDROOT/.git" ]; then

    GITBRANCH=$(git --git-dir "$AEROSOLVEDROOT/.git" rev-parse --abbrev-ref HEAD 2>/dev/null)
    GITCOMMIT=$(git --git-dir "$AEROSOLVEDROOT/.git" rev-parse HEAD 2>/dev/null)

    # Capture the current diff-stat (empty when tree is clean).
    git --git-dir "$AEROSOLVEDROOT/.git" diff --stat=999 > "$TEMPGITINFOFILE" 2>/dev/null

    # Composite signature: branch | commit | diffstat (so uncommitted changes
    # in the tree also trigger a regenerate, matching upstream behaviour).
    SIGNATURE=$(printf '%s|%s|%s' "$GITBRANCH" "$GITCOMMIT" "$(cat "$TEMPGITINFOFILE")")

    if [ -f "$COMMITFILE" ] && [ "$(cat "$COMMITFILE")" = "$SIGNATURE" ]; then
        NEEDS=0
    fi

    if [ "$NEEDS" -eq 1 ]; then

        rm -f "$GITINFOFILE"

        # Emit body lines via printf '%s\n' — no backslash interpretation,
        # so the .H gets exactly the C++ source we intend.
        printf '%s\n' "void Foam::aerosolModel::gitInfo()"       >  "$GITINFOFILE"
        printf '%s\n' "{"                                          >> "$GITINFOFILE"
        printf '%s\n' '    Info<< "Git state at compilation:" << nl << nl'    >> "$GITINFOFILE"
        printf '%s\n' '         << "    Git branch = " << "'"$GITBRANCH"'" << nl'         >> "$GITINFOFILE"
        printf '%s\n' '         << "    Git commit = " << "'"$GITCOMMIT"'" << nl;'       >> "$GITINFOFILE"

        if [ -s "$TEMPGITINFOFILE" ]; then
            printf '%s\n' '    Info<< "    Git diff =" << nl;'       >> "$GITINFOFILE"
             # One complete `Info<< ... << nl;` statement per diff line — this
             # matches upstream's structure and is valid C++ (each line is a
             # self-contained expression statement; a leading-bare `<<` on its
             # own line would not compile).
            while IFS= read -r LINE; do
                 # Skip truly empty diff lines.
                if [ -n "$LINE" ]; then
                    printf '%s\n' '    Info<< "    '"$LINE"'" << nl;' >> "$GITINFOFILE"
                fi
            done < "$TEMPGITINFOFILE"
        else
            printf '%s\n' '    Info<< "    No differences found w.r.t. this commit" << nl;' >> "$GITINFOFILE"
        fi

        printf '%s\n' "}"  >> "$GITINFOFILE"

        # Record the signature so repeated builds skip regeneration.
        printf '%s' "$SIGNATURE" > "$COMMITFILE"

        if [ -f "$TEMPGITINFOFILE" ]; then
            rm -f "$TEMPGITINFOFILE"
        fi
    fi

else

    # git absent or not a repo → single static stub.  Idempotent against
    # repeated invocations via the "no-git" marker sidecar.
    if [ ! -f "$COMMITFILE" ] || [ "$(cat "$COMMITFILE")" != "no-git" ]; then
        rm -f "$GITINFOFILE"
        printf '%s\n' "void Foam::aerosolModel::gitInfo()"   >  "$GITINFOFILE"
        printf '%s\n' "{"                                      >> "$GITINFOFILE"
        printf '%s\n' '    Info<< "    Git is not installed or this is not a git repository" << nl;' >> "$GITINFOFILE"
        printf '%s\n' "}"                                       >> "$GITINFOFILE"
        printf '%s' "no-git" > "$COMMITFILE"
    fi
fi

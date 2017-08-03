#!/usr/bin/env bash

_shogun_completion() {
    COMPREPLY=( $( env COMP_WORDS="${COMP_WORDS[*]}" \
                   COMP_CWORD=$COMP_CWORD \
                   _SHOGUN_COMPLETE=complete $1 ) )
    return 0
}

complete -F _shogun_completion -o default shogun;

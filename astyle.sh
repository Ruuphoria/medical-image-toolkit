
#!/bin/bash
astyle \
        --style=allman \
        \
        --indent=spaces=4  \
        --indent-switches \
        --indent-col1-comments \
        \
        --pad-oper \
        --pad-paren-out \
        --pad-header \
        \
        --convert-tabs \
        --break-closing-brackets \
        --add-brackets \
        --align-pointer=type \
        --max-instatement-indent=80 \
        --preserve-date \
        --lineend=linux \
        \
        --suffix=none \
        --verbose \
        --recursive "*.*xx"
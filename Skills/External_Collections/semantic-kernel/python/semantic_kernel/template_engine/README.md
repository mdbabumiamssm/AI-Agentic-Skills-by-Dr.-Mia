<!--
# COPYRIGHT NOTICE
# This file is part of the "Universal Biomedical Skills" project.
# Copyright (c) 2026 MD BABU MIA, PhD <md.babu.mia@mssm.edu>
# All Rights Reserved.
#
# This code is proprietary and confidential.
# Unauthorized copying of this file, via any medium is strictly prohibited.
#
# Provenance: Authenticated by MD BABU MIA

-->

# Prompt Template Engine

The Semantic Kernel uses the following grammar to parse prompt templates:

```
# BNF parsed by TemplateTokenizer
[template]       ::= "" | [block] | [block] [template]
[block]          ::= [sk-block] | [text-block]
[sk-block]       ::= "{{" [variable] "}}" | "{{" [value] "}}" | "{{" [function-call] "}}"
[text-block]     ::= [any-char] | [any-char] [text-block]
[any-char]       ::= any char

# BNF parsed by CodeTokenizer:
[template]       ::= "" | [variable] " " [template] | [value] " " [template] | [function-call] " " [template]
[variable]       ::= "$" [valid-name]
[value]          ::= "'" [text] "'" | '"' [text] '"'
[function-call]  ::= [function-id] | [function-id] [parameter]
[parameter]      ::= [variable] | [value]

# BNF parsed by dedicated blocks
[function-id]    ::= [valid-name] | [valid-name] "." [valid-name]
[valid-name]     ::= [valid-symbol] | [valid-symbol] [valid-name]
[valid-symbol]   ::= [letter] | [digit] | "_"
[letter]         ::= "a" | "b" ... | "z" | "A" | "B" ... | "Z"
[digit]          ::= "0" | "1" | "2" | "3" | "4" | "5" | "6" | "7" | "8" | "9"
```


<!-- AUTHOR_SIGNATURE: 9a7f3c2e-MD-BABU-MIA-2026-MSSM-SECURE -->
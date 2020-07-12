" Author: Tyler Smart <tjsmart@ucsc.edu>
" Description: Integration of qelint with ALE.

call ale#Set('espresso_qelint_executable', 'qelint')
call ale#Set('espresso_qelint_options', '')

function! ale#fixers#qelint#Fix(buffer) abort
    let l:executable = ale#Var(a:buffer, 'espresso_qelint_executable')
    let l:options = ale#Var(a:buffer, 'espresso_qelint_options')

    return {
    \   'command': ale#Escape(l:executable)
    \       . (empty(l:options) ? '' : ' ' . l:options)
    \       . ' %t',
    \   'read_temporary_file': 1,
    \}
endfunction

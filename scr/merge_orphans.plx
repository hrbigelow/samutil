#!/usr/bin/perl



@FH = [$F1, $F2];

for $f (1..2)
{
    open F, '<' . $ARGV[$f-1];
    while ($id = <F>)
    {
        $seq = <F>;
        <F>;
        $qual = <F>;
        chomp $id;
        chomp $seq;
        chomp $qual;
        $S{$id}{$f} = [$seq, $qual];
    }
    close F;
}


open OF1, ">$ARGV[2]";
open OF2, ">$ARGV[3]";

for $id (keys %S)
{
    if (scalar keys %{$S{$id}} == 2)
    {
        printf OF1 "%s\n%s\n+\n%s\n", $id, $S{$id}{1}[0], $S{$id}{1}[1];
        printf OF2 "%s\n%s\n+\n%s\n", $id, $S{$id}{2}[0], $S{$id}{2}[1];
    }
}

close OF1;
close OF2;

use ExtUtils::MakeMaker;

WriteMakefile(
    NAME              => "PrimerMapper",
    VERSION_FROM      => "lib/PrimerMapper.pm", # finds $VERSION
    PREREQ_PM         => {
        "Tkx"                         => 0,
        "LWP::Simple"                 => 0,
        "LWP::UserAgent"              => 0,
        "URI::Escape"                 => 0,
        "HTTP::Request::Common"       => 0,
        "Bio::SeqIO"                  => 0,
        "Bio::Graphics"               => 0,
        "Bio::SeqFeature::Generic"    => 0,
    },
    TEST_REQUIRES     => {
        "Test::More"                  => 0,
    },
);

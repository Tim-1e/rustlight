thanks for your reply,but something weird happened to the reference,as shown in the picture below compared with ualpha.png and smis_all_k2_stratified.png that I render by the example command line:

it seems the lightness turned to be darker

As for the comparison format as png,because I meet save problem when I use exr format,
the code `img.save("primal", imgout_path_str);`causes errors for exr format,and I don't know exactly it's the lib related problem due to the environment or rust's version problem,which I might figure out later.

And about the "pbrt render problem",

> Additionally,I am experiencing difficulties in re-rendering the PBRT scene using PBRT-v3. The rendered images do not seem to match the expected results.So I either can get the correct result by this way.

I use the scene you mentioned in the example as [scene](https://data.adrien-gruson.com/research/2020_CMIS/plane_scene.zip),which has pbrt file.
I use the `meeting_ply.pbrt` in the unzip filefold for pbrt-v3,but I found it didn't match the input format of pbrt-v3 as it lack of `WorldBegin` sign,
so I add a `WorldBegin` sign before the `AttributeBegin` sign and it could run through the pbrt-v3,but the result is weird as shown in the picture below:


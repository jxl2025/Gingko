#include "sceneparser.h"
#include "scenefilereader.h"
#include <glm/gtx/transform.hpp>

#include <chrono>
#include <iostream>

struct Word {
    std::string word;
    Word* left;
    Word* right;
};

void dfsPrintTree(Word* w, std::string sentence) {
    // Task 4: Debug this function!!! (Hint: you may need to write/add some logic...)
    // std::string newSentence = w->word + sentence;
    // if (w->left == nullptr && w->right == nullptr) {
    //     std::cout << newSentence << std::endl;
    // }
    // dfsPrintTree(w->left, newSentence);
    // dfsPrintTree(w->right, newSentence);
    // two problems in above code:
    // 1. if we dereference null pointer, we get error
    // 2. we need to stop at leaves (rather than keep recursing

    if (!w) return;
    std::string newSentence;
    if (sentence.empty()) newSentence = w->word;
    else newSentence = sentence + " " + w->word;

    // Leaf: no children → print
    if (!w->left && !w->right) {
        std::cout << newSentence << std::endl;
        return;
    }

    // Recurse down both sides
    dfsPrintTree(w->left,  newSentence);
    dfsPrintTree(w->right, newSentence);
}

Word* initTree(std::vector<Word> &words) {
    // STUDENTS - DO NOT TOUCH THIS FUNCTION
    words.reserve(8);
    words.push_back(Word{"2D graphics ", nullptr, nullptr});
    words.push_back(Word{"3D graphics ", nullptr, nullptr});
    words.push_back(Word{"making ", &words[0], &words[1]});
    words.push_back(Word{"CS123 ", nullptr, nullptr});
    words.push_back(Word{"love ", &words[3], &words[2]});
    words.push_back(Word{"bugs ", nullptr, nullptr});
    words.push_back(Word{"hate ", nullptr, &words[5]});
    words.push_back(Word{"I ", &words[4], &words[6]});
    return &words[7];
}


void traverseNode(SceneNode* node, const glm::mat4& parentCTM,
                  std::vector<RenderShapeData>& shapes, std::vector<SceneLightData>& lights)
{
    if (!node) return;
    // local transform in order listed
    glm::mat4 M = parentCTM;
    for (const SceneTransformation* t : node->transformations) {
        switch (t->type) {
        case TransformationType::TRANSFORMATION_TRANSLATE:
            M *= glm::translate(glm::mat4(1.f), t->translate);
            break;
        case TransformationType::TRANSFORMATION_SCALE:
            M *= glm::scale(glm::mat4(1.f), t->scale);
            break;
        case TransformationType::TRANSFORMATION_ROTATE:
            // angle is already in radians per ScenefileReader
            M *= glm::rotate(glm::mat4(1.f), t->angle, t->rotate);
            break;
        case TransformationType::TRANSFORMATION_MATRIX:
            M *= t->matrix;
            break;
        }
    }
    // get the primitives
    for (const ScenePrimitive* p : node->primitives) {
        shapes.push_back(RenderShapeData{ *p, M });}
    //lighting
    for (const SceneLight* L : node->lights) {
        SceneLightData W{};
        W.id = L->id;
        W.type = L->type;
        W.color = L->color;
        W.function  = L->function;
        W.penumbra = L->penumbra;
        W.angle = L->angle;
        W.width = L->width;
        W.height= L->height;
        // transform node origin (need for directional, I think?)
        if (W.type != LightType::LIGHT_DIRECTIONAL) {
            W.pos = M * glm::vec4(0.f, 0.f, 0.f, 1.f);
        }
        // transform direction and normalize it
        if (W.type != LightType::LIGHT_POINT) {
            glm::vec4 dirWorld = M * glm::vec4(glm::vec3(L->dir), 0.f);
            W.dir = glm::vec4(glm::normalize(glm::vec3(dirWorld)), 0.f);
        }
        lights.push_back(W);
    }
    for (SceneNode* c : node->children) {
        traverseNode(c, M, shapes, lights);
    }
}


bool SceneParser::parse(std::string filepath, RenderData &renderData) {
    ScenefileReader fileReader = ScenefileReader(filepath);
    bool success = fileReader.readJSON();
    if (!success) {
        return false;
    }

    // TODO: Use your Lab 5 code here
    // Task 5: populate renderData with global data, and camera data;
    renderData.globalData = fileReader.getGlobalData();
    renderData.cameraData = fileReader.getCameraData();
    // Task 6: populate renderData's list of primitives and their transforms.
    //         This will involve traversing the scene graph, and we recommend you
    //         create a helper function to do so!
    renderData.shapes.clear();
    renderData.lights.clear();
    SceneNode* root = fileReader.getRootNode();
    traverseNode(root, glm::mat4(1.f), renderData.shapes, renderData.lights);

    return true;
}

